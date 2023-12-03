import numpy as np

def date_finder(j):
        if j < 31:
            mon = "1912"
            dayint = j + 1
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)
        if j > 30 and j < 62:
            mon = "2001"
            dayint = j - 30
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)
        if j > 61 and j < 91:
            mon = "2002"
            dayint = j - 61
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)
        if j > 90 and j < 122:
            mon = "2003"
            dayint = j - 90
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)
        if j > 121 and j < 152:
            mon = "2004"
            dayint = j - 121
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)
        if j > 151 and j < 183:
            mon = "2005"
            dayint = j - 151
            if dayint < 10:
                day = "0"+str(dayint)
            if dayint > 9:
                day = str(dayint)

        return day, mon
