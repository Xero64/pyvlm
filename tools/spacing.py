from math import pi, cos, asin

def normalise_spacing(spacing: list) -> list:
    s = sorted(spacing)
    smin = min(s)
    smax = max(s)
    news = [(sval-smin)/(smax-smin) for sval in s]
    return news

def semi_elliptical_spacing(num: int) -> list:
    s = [cos(asin(1-float(i)/num)) for i in range(num+1)]
    s[0] = 0.0
    return s

def full_elliptical_spacing(num: int) -> list:
    th = [pi-asin(float(i)/num*2) if i < num/2 else asin(2*(num-float(i))/num) for i in range(num+1)]
    s = [(1+cos(thi))/2 for thi in th]
    return s

def semi_cosine_spacing(num: int) -> list:
    s = [cos(float(num-i)*pi/2/num) for i in range(num+1)]
    s[0] = 0.0
    return s

def full_cosine_spacing(num: int) -> list:
    s = [(1+cos(float(num-i)*pi/num))/2 for i in range(num+1)]
    return s

def equal_spacing(num: int) -> list:
    return [float(i)/num for i in range(num+1)]
