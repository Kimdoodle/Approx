#include "header/CKKS_params.h"

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

CKKS_params::CKKS_params(vector<int> modulus, double scale, size_t pmd, size_t hwt)
{
    this->scale = scale;
    parms = make_unique<EncryptionParameters>(scheme_type::ckks);
    parms->set_poly_modulus_degree(pmd);
   
    parms->set_coeff_modulus(CoeffModulus::Create(pmd, modulus));
    
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
// encode coeff, but scale, parms_id will be equal to ctxt's.
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

// encode+encrypt plain value, but params will be equal with ctxt's.
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

// Add
void CKKS_params::add(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    scale_equal(ctxt1, ctxt2);
    scale_equal(ctxt1, result);
    eva->add(ctxt1, ctxt2, result);
}
void CKKS_params::add(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    scale_equal(ctxt1, ctxt2);
    eva->add_inplace(ctxt1, ctxt2);
}
void CKKS_params::add(Plaintext& ptxt, Ciphertext& ctxt)
{
    //scale_equal(ptxt, ctxt);
    eva->add_plain_inplace(ctxt, ptxt);
}

// Multiply
void CKKS_params::mult(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result)
{
    //modulus_equal(ctxt1, ctxt2);
    scale_equal(ctxt1, ctxt2);
    eva->multiply(ctxt1, ctxt2, result);
    eva->relinearize_inplace(result, rlk);
    eva->rescale_to_next_inplace(result);
}

//Multiply 2 ciphertexts.
void CKKS_params::mult(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    scale_equal(ctxt1, ctxt2);
    eva->multiply_inplace(ctxt1, ctxt2);
    eva->relinearize_inplace(ctxt1, rlk);
    eva->rescale_to_next_inplace(ctxt1);
}

//Multiply plaintext / ciphertext. two data's scale, modulus_level should be equal.
void CKKS_params::mult(Plaintext& ptxt, Ciphertext& ctxt)
{
    eva->multiply_plain_inplace(ctxt, ptxt);
    eva->rescale_to_next_inplace(ctxt);
}

//square ciphertext.
void CKKS_params::square(Ciphertext& ctxt)
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
            mult(result, squareX);
        }
        square(squareX);
        d /= 2;
    }
    return result;
}

//두 암호문의 scale 통일
void CKKS_params::scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2)
{
    while (ctxt1.coeff_modulus_size() > ctxt2.coeff_modulus_size()) {
        Plaintext p = encode(1.0, ctxt1);
        mult(p, ctxt1);
    }

    while (ctxt1.coeff_modulus_size() < ctxt2.coeff_modulus_size()) {
        Plaintext p = encode(1.0, ctxt2);
        mult(p, ctxt2);
    }
}
void CKKS_params::scale_equal(Plaintext& ptxt, Ciphertext& ctxt)
{
    Ciphertext x;
    enc->encrypt(ptxt, x);
    scale_equal(x, ctxt);
    dec->decrypt(x, ptxt);
}
