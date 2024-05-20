#####------------------------ Useful Functions -----------------------#####
def subcirc(l1, l2):
    return [(l1[i]-l2[i]+180)%360-180 for i in range(len(l1))]

def subcirc_mins(l1, l2):
    return [60*((l1[i]-l2[i]+180)%360-180) for i in range(len(l1))]

def subcirc_secs(l1, l2):
    return [3600*((l1[i]-l2[i]+180)%360-180) for i in range(len(l1))]

def todms(_angle):
    if _angle > 0:
        sign = 1
    else:
        _angle = -_angle
        sign = -1
    ras = _angle/30
    deg = 30*(ras - int(ras))
    min = 60*(deg - int(deg))
    sec = 60*(min - int(min))
    thd = 60*(sec - int(sec))
    return [sign*int(ras), int(deg), int(min), int(sec), thd]

def limsmall(angle):
    return (angle+180)%360-180

def intToRoman(num):
    # Storing roman values of digits from 0-9
    # when placed at different places
    m = ["", "M", "MM", "MMM"]
    c = ["", "C", "CC", "CCC", "CD", "D", "DC", "DCC", "DCCC", "CM"]
    x = ["", "X", "XX", "XXX", "XL", "L", "LX", "LXX", "LXXX", "XC"]
    i = ["", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX"]
    # Converting to roman
    thousands = m[num // 1000]
    hundreds = c[(num % 1000) // 100]
    tens = x[(num % 100) // 10]
    ones = i[num % 10]
 
    ans = (thousands + hundreds + tens + ones) 
    return ans