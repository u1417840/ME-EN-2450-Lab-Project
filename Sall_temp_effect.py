def sall_temp_effect(T):
    if T <= 0:
        PT = 0
    elif 0 < T < 35:
        PT = 0.000241 * T**2.06737 * (35 - T)**0.72859
    else:
        PT = 0
    return PT
