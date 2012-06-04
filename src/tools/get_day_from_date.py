"""
returns numbers of days counted from 01JAN2001 as day 0

date must be provided in format
DDMMMYYY where MMM is taken from
["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"]

negative values are probabely wrong
"""

def get_day_from_date(dstr):
    monthdict = {"JAN":1,"FEB":2,"MAR":3,
                 "APR":4,"MAY":5,"JUN":6,
                 "JUL":7,"AUG":8,"SEP":9,
                 "OCT":10,"NOV":11,"DEC":12}

    list_of_days_per_month = [31,28,31,30,31,30,31,31,30,31,30,31]
    days_per_month = dict([(i+1, sum(list_of_days_per_month[0:i])) for i in range(12)])

    day = int(dstr[0:2])
    month = monthdict[dstr[2:5]]
    year = int(dstr[5:9])

    days = (year - 2001) * 365
    days += days_per_month[month]
    days += day
    #schaltjahre
    days += (year - 2005)/4
    if (year - 2004)%4 == 0 and month >= 3:
        days += 1
    days -= 1 #index starts at zero
    return days

if __name__ == "__main__":
    print get_day_from_date("01JAN2008")

