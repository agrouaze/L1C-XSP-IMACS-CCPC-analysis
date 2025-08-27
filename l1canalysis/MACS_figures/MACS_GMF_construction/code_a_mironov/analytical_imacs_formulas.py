import numpy as np


# y = (((exp(cos3(x₂)) + ((x₁ * 1.0981) - 6.099)) * (cos3(x₀ ^ -1.9469) - 1.237)) * (cos(x₂ + -0.0175) + 0.052687)) + (tan(sin(cos2(x₂))) * sin((x₁ * x₀) * 0.71937))
def imacs_fit9_agrouaze(alpha, u, phi):
    """

    :param alpha : (np.ndarray) incidences angles
    :param u: (np.ndarray) wind speed
    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    # 10 jan 2025
    alpha = np.deg2rad(alpha)
    fit = (
        (
            (np.exp(np.cos(phi) ** 3) + ((u * 1.0981) - 6.099))
            * (np.cos(alpha ** -1.9469) ** 3 - 1.237)
        )
        * (np.cos(phi + -0.0175) + 0.052687)
    ) + (np.tan(np.sin(np.cos(phi) ** 2)) * np.sin((u * alpha) * 0.71937))
    imacs = fit / 100
    return imacs

#35         loss= 6.332e-02  score=2.017e-03  y = (((cos2(x₂ - 0.016237) * (sin(exp(x₀) * (x₁ / 4.241)) + -1.1348)) + (cos(x₂ + -0.020324) * (4.3902 - x₁))) * cos3(sin(cos3(tan(-1.3652 - x₀)) - -0.62295))) * 2.3921 
#/raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Feb11_0558.pkl # run with 40k iterations before the change of azimuths convention and with the simplification of gridata+RectBivariate
# I think the x variable was wdir_az_scat (0 - 360)
def imacs_fit10_agrouaze(phi,alpha, u,c1=0.016237,c2=4.241,c3=-1.1348,c4=-0.020324,c5=4.3902,c6=-1.3652,c7=-0.62295,c8=2.3921):
    """

    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise
    :param alpha : (np.ndarray) incidences angles
    :param u: (np.ndarray) wind speed
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    # 10 jan 2025
    alpha = np.deg2rad(alpha)
    fit = (((np.cos(phi - c1)**2 * (np.sin(np.exp(alpha) * (u / c2)) + c3)) + (np.cos(phi + c4) * (c5 - u))) * np.cos(np.sin(np.cos(np.tan(c6 - alpha))**3 - c7))**3) * c8
    imacs = fit / 100
    return imacs

# 28          7.303e-02  8.617e-02  y = (((sin(x₁ / exp(cos(x₀))) + -1.1079) * cos2(x₂)) + (cos(x₂) * (4.4047 - x₁))) * (1.4088 - cos3(tan(x₀ - -1.2719) - -0.45106))
# /raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Feb11_0558.pkl # run with 40k iterations before the change of azimuths convention and with the simplification of gridata+RectBivariate
def imacs_fit11_agrouaze(phi,alpha, u,c1=-1.1079,c2=4.4047,c3=1.4088,c4=-1.2719,c5=-0.45106):
    """

    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2
    :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
    :param u: (np.ndarray) wind speed x1
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    # 10 jan 2025
    alpha = np.deg2rad(alpha)
    fit = (((np.sin(u / np.exp(np.cos(alpha))) + c1) * np.cos(phi)**2) + (np.cos(phi) * (c2 - u))) * (c3 - np.cos(np.tan(alpha - c4) - c5)**3)
    imacs = fit / 100
    return imacs

# 35          1.799e-01  1.231e-02  y = (inv(x₀ - 0.34032) * ((((4.773 - x₁) * tan(cos(x₂ + -0.032276))) * 0.39568) - cos2(tan(cos2(cos(x₂) - -0.23528))))) - ((cos2(x₂) * cos(sin(x₀) * x₁)) + -0.3432)
# /raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Mar04_1439.pkl # run with 40k iterations after change of azimuths convention and with gaussian filter 1D
def imacs_fit12_agrouaze(phi,alpha, u):
    """

        :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2
        :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
        :param u: (np.ndarray) wind speed x1
        :return:
            imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    # 6 march 2025
    alpha = np.deg2rad(alpha)
    fit = (1/(alpha - 0.34032) * ((((4.773 - u) * np.tan(np.cos(phi + -0.032276))) * 0.39568) - np.cos(
        np.tan(np.cos(np.cos(phi) - -0.23528)**2))**2)) - ((np.cos(phi)**2 * np.cos(np.sin(alpha) * u)) + -0.3432)
    imacs = fit / 100
    return imacs

# 35          5.503e-01  5.217e-02  y = ((((cos(x₂) + (x₁ + -5.3842)) / (0.31512 - x₀)) - cos(x₁)) + cos(tan(x₀ - 2.0622))) * (tan(cos(x₂) - ((cos(x₁ * sin(x₀)) * -0.047966) + 0.041049)) - -0.11887)
# /raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Mar12_1002.pkl training IMACS pyssr méthode gaussian filter (plage de vent limité) lambda max = 100m
def imacs_fit13_agrouaze(phi,alpha,u):
    """

            :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
            :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
            :param u: (np.ndarray) wind speed x1
            :return:
                imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

        """
    # 12 march 2025
    alpha = np.deg2rad(alpha)
    fit = ((((np.cos(phi) + (u + -5.3842)) / (0.31512 - alpha)) - np.cos(u)) + np.cos(np.tan(alpha - 2.0622))) * (np.tan(np.cos(phi) - ((np.cos(u * np.sin(alpha)) * -0.047966) + 0.041049)) - -0.11887)
    imacs = fit / 100
    return imacs

#35          5.466e-01  1.588e-02  y = tan(cos(x₂ - 0.021366)) * (((sin(exp(x₀) * (cos(x₂) - (x₁ / 6.3776))) - (sin(x₀) * (x₁ + -5.5862))) / (-0.43091 + x₀)) - ((cos(x₁) - -1.8022) ^ cos2(x₂)))
# /raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Mar22_0941.pkl  training IMACS pyssr méthode gaussian filter (plage de vent limité) lambda max = 100m, variable sigma
def imacs_fit14_agrouaze(phi,alpha,u):
    """

            :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
            :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
            :param u: (np.ndarray) wind speed x1
            :return:
                imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

        """
    # 22 march 2025
    alpha = np.deg2rad(alpha)
    fit =  np.tan(np.cos(phi - 0.021366)) * (((np.sin(np.exp(alpha) * (np.cos(phi) - (u / 6.3776))) - (np.sin(alpha) * (u + -5.5862))) / (-0.43091 + alpha)) - ((np.cos(u) - -1.8022) ** np.cos(phi)**2))
    imacs = fit / 100
    return imacs


#34 complexity ,(-1.6030896*cos(x2 - 0.042280238)**2 - 0.30620417/(x0 - 0.46252897))*((x1 - 3.145058)*(cos(x2) - 0.10230599*tan(x0)) - 1.9961689*sin(0.4329867*x1 - cos(x2))*cos(x2)**2 + 1.3950036)
#/raid/localscratch/agrouaze/imacs_gmf/iw/pysr_fit_run_2025Apr07_1053.pkl 1000x20 iteration with imacs_nn_fit_adapted_from_mironov_grougrou3.py, loss stuck at value 5 for val and train.
def imacs_fit15_agrouaze(phi,alpha,u):
    """

            :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
            :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
            :param u: (np.ndarray) wind speed x1
            :return:
                imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

        """
    # 22 march 2025
    alpha = np.deg2rad(alpha)
    fit = (-1.6030896*np.cos(phi - 0.042280238)**2 - 0.30620417/(alpha - 0.46252897))*((u - 3.145058)*(np.cos(phi) - 0.10230599*np.tan(alpha)) - 1.9961689*np.sin(0.4329867*u - np.cos(phi))*np.cos(phi)**2 + 1.3950036)
    imacs = fit / 100
    return imacs

# first training on residual IMACS based on Acos(phi)+residual, 14 May 2025
# http://134.246.184.23:6080/#/experiments/878239617694006373/runs/7186b8a0989749768d36fe2559989879


def imacsAF16(phi,alpha,u):
    """
    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
    :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
    :param u: (np.ndarray) wind speed x1
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    alpha = np.deg2rad(alpha)
    # residual equ/ y = ((((4.2422 - x₁) / x₀) * sin(x₂ - 1.5775)) + 0.89615) * ((-0.98418 / x₀) - cos2(x₂))
    # fit = (np.cos(np.radians(phi)) * 3.1333) * (4.0298 - u) +  ((((4.2422 - u) / alpha) * np.sin(phi - 1.5775)) + 0.89615) * ((-0.98418 / alpha) - np.cos(phi)**2)
    fit = (np.cos(phi) * 3.1333) * (4.0298 - u) + (
                (((4.2422 - u) / alpha) * np.sin(phi - 1.5775)) + 0.89615) * ((-0.98418 / alpha) - np.cos(phi) ** 2)
    imacs = fit / 100
    return imacs

# y = ((-1.0541 - (tan(cos(x₂)) * x₁)) / tan(x₀ + cos2((exp(cos(x₂)) + x₁) / 10.593))) / x₀
# best equation
# /raid/localscratch/agrouaze/imacs_gmf/iw/2025May16-1541/20250516_154116_xz1oEP/hall_of_fame.csv
# 21 may 2025 fun-pug-156 http://134.246.184.23:6080/#/experiments/878239617694006373/runs/e9229fa8731b4e5f8664b72a9f66c454
def imacsAF17(phi,alpha,u):
    """

    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
    :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
    :param u: (np.ndarray) wind speed x1
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    alpha = np.deg2rad(alpha)
    # residual equ/ y = ((((4.2422 - x₁) / x₀) * sin(x₂ - 1.5775)) + 0.89615) * ((-0.98418 / x₀) - cos2(x₂))
    epsilon = ((-1.0541 - (np.tan(np.cos(phi)) * u)) / np.tan(alpha + np.cos((np.exp(np.cos(phi)) + u) / 10.593)**2)) / alpha
    fit = (np.cos(phi) * 3.1333) * (4.0298 - u) + epsilon/100
    print('sans epsilon')
    imacs = fit / 100
    return imacs

# y = ((cos(x₂) * ((x₁ + -4.2225) - sin(x₁ * -0.56163))) - tan(cos((x₁ / (-2.6992 / x₀)) + -0.72456))) * ((6.5571 - ((x₀ * 1.2769) ^ x₁)) + (cos(x₂ + (x₂ + 3.0573)) + ...
#                                       (-6.061 / x₀)))
# longest equation (higher complexity)
# /raid/localscratch/agrouaze/imacs_gmf/iw/2025May16-1541/20250516_154116_xz1oEP/hall_of_fame.csv
# 21 may 2025 fun-pug-156 http://134.246.184.23:6080/#/experiments/878239617694006373/runs/e9229fa8731b4e5f8664b72a9f66c454
def imacsAF18(phi,alpha,u):
    """

    :param phi: (np.ndarray) azimuth wind direction (relative to antenna), clockwise degrees x2 [radians]
    :param alpha : (np.ndarray) incidences angles x0 (in the pyssr training order)
    :param u: (np.ndarray) wind speed x1
    :return:
        imacs (np.ndarray): IMACS VV intraburst GMF-like 17.5km² trained from NN analytical formula

    """
    alpha = np.deg2rad(alpha)
    # residual equ/ y = ((((4.2422 - x₁) / x₀) * sin(x₂ - 1.5775)) + 0.89615) * ((-0.98418 / x₀) - cos2(x₂))
    epsilon = ((np.cos(phi) * ((u + -4.2225) - np.sin(u * -0.56163))) - np.tan(np.cos((u / (-2.6992 / alpha)) + -0.72456))) * ((6.5571 - ((alpha * 1.2769) ** u)) + (np.cos(phi + (phi + 3.0573)) + (-6.061 / alpha)))
    fit = (np.cos(phi) * 3.1333) * (4.0298 - u) + epsilon
    imacs = fit / 100
    return imacs