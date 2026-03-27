/*
    현재 상황
    Decomp -> EvalStep : python에서 수행
    EvalStep -> Boundary : python에서 c++코드로 수행
    
    목표
    coeff가 제시될 때 
    Decomp(cache)불러오기 -> EvalStep -> Boundary
    모두 python에서 실행하지만 c++코드로 수행

    1. Decomp(cache)
    - json.hpp를 사용하여 load_decomp_cache를 호출.
    - reconstruct_decomp_from_cache 호출.
    
    2. EvalStep
    - Decomp만 잘 생성하면 생성자로 간단하게 변환가능.

    3. Boundary
    - fct::Ciphertext::evaluate로 수행가능.
*/