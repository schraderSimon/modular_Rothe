def get_individual_params(params):
    a, b, mu, p = params.transpose(2, 0, 1)
    return a, b, mu, p
