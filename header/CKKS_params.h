#pragma once
#include "seal/seal.h"
#include "seal/bootstrapper.h"
using namespace std;
using namespace seal;


class CKKS_params {
public:
    unique_ptr<EncryptionParameters> parms;
    vector<int> modulus;
    unique_ptr<SEALContext> context;
    unique_ptr<CKKSBootstrapper> boot;
    unique_ptr<KeyGenerator> keygen;
    SecretKey sk;
    PublicKey pk;
    RelinKeys rlk;

    unique_ptr<Encryptor> enc;
    unique_ptr<Evaluator> eva;
    unique_ptr<Decryptor> dec;
    unique_ptr<CKKSEncoder> encoder;

    double scale;

    CKKS_params(vector<int> modulus, double scale, size_t pmd, size_t hwt);

    Plaintext encode(double input);
    Plaintext encode(const double& input, Ciphertext& ctxt);
    Plaintext encode(const double& input, Ciphertext& ctxt, double scale);
    Plaintext encode(const vector<double>& input);
    Ciphertext encrypt(const Plaintext& plain);
    Ciphertext encrypt(double input);
    Ciphertext encrypt(const double& input, Ciphertext& ctxt);
    Ciphertext encrypt(const vector<double>& input);

    Plaintext decrypt(Ciphertext& ctxt);
    vector<double> decode(const Plaintext& ptxt);
    vector<double> decode_ctxt(Ciphertext& ctxt);

    void modulus_switch(Plaintext& ptxt, const parms_id_type parms_id);
    void modulus_switch(Ciphertext& ctxt, const parms_id_type parms_id);
    void add(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result);
    void add(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void add(Plaintext& ptxt, Ciphertext& ctxt);
    void mult(Ciphertext& ctxt1, Ciphertext& ctxt2, Ciphertext& result);
    void mult(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void mult(Plaintext& ptxt, Ciphertext& ctxt);
    void square(Ciphertext& ctxt);

    Ciphertext exp(const Ciphertext& x, int d);
    void modulus_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Ciphertext& ctxt1, Ciphertext& ctxt2);
    void scale_equal(Plaintext& ptxt, Ciphertext& ctxt);
};