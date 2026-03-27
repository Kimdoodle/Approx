#include "../include/CKKS_params.h"

void print_parameters(const SEALContext& context)
{
    auto& context_data = *context.key_context_data();
    string scheme_name = "CKKS";

    cout << "| Encryption parameters :" << endl;
    cout << "|   scheme: " << scheme_name << endl;
    cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    cout << "|   coeff_modulus size: ";
    cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    size_t coeff_modulus_size = coeff_modulus.size();
    for (size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        cout << coeff_modulus[i].bit_count() << " + ";
    }
    cout << coeff_modulus.back().bit_count();
    cout << ") bits" << endl;
    cout << "---" << endl;
}

CKKS_params::CKKS_params(vector<int> modulus, int p_num, int e_num, int s_num, int N, size_t hwt, pair<int, int> iter_HELUT)
{
    this->p = p_num, this->e = e_num, this->s = s_num, this->N = N;
    this->x_bound = ldexp(1.0, p_num);
    this->target_precision = ldexp(1.0, -e_num);
    this->scale = ldexp(1.0, s_num);
    this->poly_modulus_degree = size_t(1) << N;

    this->iter_HELUT = iter_HELUT;
    parms = make_unique<EncryptionParameters>(scheme_type::ckks);
    parms->set_poly_modulus_degree(this->poly_modulus_degree);
   
    parms->set_coeff_modulus(CoeffModulus::Create(this->poly_modulus_degree, modulus));
    
    context = make_unique<SEALContext>(*parms, true, sec_level_type::none);
    //boot = make_unique<CKKSBootstrapper>(*context);
    //print_parameters(*context);
    keygen = make_unique<KeyGenerator>(*context, hwt);
    sk = keygen->secret_key();
    keygen->create_public_key(pk);
    keygen->create_relin_keys(rlk);

    enc = make_unique<Encryptor>(*context, pk);
    eva = make_unique<Evaluator>(*context);
    dec = make_unique<Decryptor>(*context, sk);
    encoder = make_unique<CKKSEncoder>(*context);
}

// encode coeff. Only 1 value needed.
Plaintext CKKS_params::encode(double input)
{
    Plaintext plain;
    encoder->encode(input, scale, plain);
    return plain;
}
// encode coeff. Scale, parms_id will be equal to ctxt's.
Plaintext CKKS_params::encode(const double& input, Ciphertext& ctxt)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), ctxt.scale(), plain);
    return plain;
}
//encode coeff, but only parms_id will be equal to ctxt's.
Plaintext CKKS_params::encode(const double& input, Ciphertext& ctxt, double scale)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), scale, plain);
    return plain;
}

// encode input value. vector needed.
Plaintext CKKS_params::encode(const vector<double>& input)
{
    Plaintext plain;
    encoder->encode(input, scale, plain);
    return plain;
}

Plaintext CKKS_params::encode(const vector<double>& input, Ciphertext& ctxt)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), ctxt.scale(), plain);
    return plain;
}

Plaintext CKKS_params::encode(const vector<double>& input, Ciphertext& ctxt, double scale)
{
    Plaintext plain;
    encoder->encode(input, ctxt.parms_id(), scale, plain);
    return plain;
}

// encrypt plaintext.
Ciphertext CKKS_params::encrypt(const Plaintext& plain)
{
    Ciphertext ctxt;
    enc->encrypt(plain, ctxt);
    return ctxt;
}

// encode+encrypt plain value.
Ciphertext CKKS_params::encrypt(double input)
{
    return encrypt(encode(input));
}

// encode+encrypt plain value. parms_id, scale will be equal to ctxt's.
Ciphertext CKKS_params::encrypt(const double& input, Ciphertext& ctxt)
{
    return encrypt(encode(input, ctxt));
}

// encode+encrypt plain vector
Ciphertext CKKS_params::encrypt(const vector<double>& input)
{
    return encrypt(encode(input));
}

// decrypt
Plaintext CKKS_params::decrypt(Ciphertext& ctxt)
{
    Plaintext dec_plain;
    dec->decrypt(ctxt, dec_plain);
    return dec_plain;
}

// decode
vector<double> CKKS_params::decode(const Plaintext& ptxt)
{
    vector<double> output;
    encoder->decode(ptxt, output);
    return output;
}

// decrypt + decode
vector<double> CKKS_params::decode_ctxt(Ciphertext& ctxt)
{
    return decode(decrypt(ctxt));
}

// Modulus Switch
void CKKS_params::modulus_switch(Plaintext& ptxt, const parms_id_type parms_id)
{
    eva->mod_switch_to_inplace(ptxt, parms_id);
}

void CKKS_params::modulus_switch(Ciphertext& ctxt, const parms_id_type parms_id)
{
    eva->mod_reduce_to_inplace(ctxt, parms_id);
}

// Addition
void CKKS_params::add_ct_ct(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    scale_equal(ctxt1, ctxt2);
    // scale_equal(ctxt1, result);
    eva->add(ctxt1, ctxt2, result);
}
void CKKS_params::add_ct_ct_inplace(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    scale_equal(ctxt1, ctxt2);
    eva->add_inplace(ctxt1, ctxt2);
}
void CKKS_params::add_pt_ct_inplace(Plaintext& ptxt, Ciphertext& ctxt)
{
    eva->add_plain_inplace(ctxt, ptxt);
}
void CKKS_params::add_pt_pt(Plaintext& pt1, Plaintext& pt2, Plaintext& result)
{
    /* Restraints
        1. both pt is ntt_form
        2. same parms_id, coeff_count, scale
    */
    if(!pt1.is_ntt_form() || !pt2.is_ntt_form())
    {
        throw std::invalid_argument("NTT form error!");
    }
    if(pt1.coeff_count() != pt2.coeff_count())
    {
        throw std::invalid_argument("coeff count error!");
    }
    if(pt1.parms_id() != pt2.parms_id())
    {
        throw std::invalid_argument("parms_id error!");
    }
    if(pt1.scale() != pt2.scale())
    {
        throw std::invalid_argument("scale error!");
    }

    auto context_data = context->get_context_data(pt1.parms_id());
    const auto &parms = context_data->parms();
    const auto &coeff_modulus = parms.coeff_modulus();
    std::size_t poly_degree = parms.poly_modulus_degree();
    std::size_t mod_count = coeff_modulus.size();
    
    result.reserve(pt1.coeff_count());
    // seal::Plaintext result(pt1.coeff_count());
    result.parms_id() = pt1.parms_id();
    result.scale() = pt1.scale();

    const auto *p1 = pt1.data();
    const auto *p2 = pt2.data();
    auto *res = result.data();
    
    for(std::size_t i=0; i<mod_count; ++i)
    {
        const auto &mod = coeff_modulus[i];
        std::size_t base = i * poly_degree;

        for(std::size_t j=0; j<poly_degree; ++j)
        {
            res[base + j] = seal::util::add_uint_mod(p1[base+j], p2[base+j], mod);
        }
    }
}
void CKKS_params::add_pt_ct(seal::Plaintext& ptxt, seal::Ciphertext& ctxt, seal::Ciphertext& res)
{
    eva->add_plain(ctxt, ptxt, res);
}

// Multiply
void CKKS_params::mult_ct_ct(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result, bool SE)
{
    // modulus_equal(ctxt1, ctxt2);
    if(SE)
        scale_equal(ctxt1, ctxt2);
    eva->multiply(ctxt1, ctxt2, result);
    eva->relinearize_inplace(result, rlk);
    eva->rescale_to_next_inplace(result);
}
void CKKS_params::mult_ct_ct_inplace(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    modulus_equal(ctxt1, ctxt2);
    // scale_equal(ctxt1, ctxt2);
    eva->multiply_inplace(ctxt1, ctxt2);
    eva->relinearize_inplace(ctxt1, rlk);
    eva->rescale_to_next_inplace(ctxt1);
}

//Multiply plaintext / ciphertext.
void CKKS_params::mult_pt_ct(seal::Plaintext& ptxt, seal::Ciphertext& ctxt, seal::Ciphertext& res, bool rescale)
{
    eva->multiply_plain(ctxt, ptxt, res);
    if(rescale)
        eva->rescale_to_next_inplace(res);
}
void CKKS_params::mult_pt_ct_inplace(Plaintext& ptxt, Ciphertext& ctxt, bool rescale)
{
    eva->multiply_plain_inplace(ctxt, ptxt);
    if(rescale)
        eva->rescale_to_next_inplace(ctxt);
}

//square ciphertext.
Ciphertext CKKS_params::square(Ciphertext& ctxt)
{
    Ciphertext result;
    eva->square(ctxt, result);
    eva->relinearize_inplace(result, rlk);
    eva->rescale_to_next_inplace(result);
    return result;
}
void CKKS_params::square_inplace(Ciphertext& ctxt)
{
    eva->square_inplace(ctxt);
    eva->relinearize_inplace(ctxt, rlk);
    eva->rescale_to_next_inplace(ctxt);
}

// exp
Ciphertext CKKS_params::exp(const Ciphertext& x, int d)
{
    Ciphertext result, squareX;
    result = encrypt(1.0);

    squareX = x;
    while(d > 0) {
        if (d % 2) {
            mult_ct_ct_inplace(result, squareX);
        }
        square(squareX);
        d /= 2;
    }
    return result;
}

//두 암호문의 modulus 통일
void CKKS_params::modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    if(ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size())
        eva->mod_reduce_to_inplace(ctxt1, ctxt2.parms_id());
    else if(ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size())
        eva->mod_reduce_to_inplace(ctxt2, ctxt1.parms_id());
}
//두 암호문의 scale 통일
void CKKS_params::scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    seal::Plaintext pt;
    while (ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size()) {
        // eva->mod_reduce_to_inplace(ctxt1, ctxt2.parms_id());
        pt = this->encode(1.0, ctxt1);
        eva->multiply_plain_inplace(ctxt1, pt);
        eva->rescale_to_next_inplace(ctxt1);
    }

    while (ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size()) {
        // eva->mod_reduce_to_inplace(ctxt2, ctxt1.parms_id());
        pt = this->encode(1.0, ctxt2);
        eva->multiply_plain_inplace(ctxt2, pt);
        eva->rescale_to_next_inplace(ctxt2);
    }
}
void CKKS_params::scale_equal(Plaintext& ptxt, Ciphertext& ctxt)
{
    Ciphertext x;
    enc->encrypt(ptxt, x);
    scale_equal(x, ctxt);
    dec->decrypt(x, ptxt);
}
