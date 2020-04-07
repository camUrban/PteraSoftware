def second_max(l):
    max_num = l[0]
    second_max_num = l[1]
    for each_num in l:
        if each_num > second_max_num and each_num <= max_num:
            second_max_num = each_num
        elif each_num >= second_max_num and each_num >= max_num:
            second_max_num = max_num
            max_num = each_num
    return second_max_num
second_max([-8, 1,  4 ,-3, 6, 8, 2, -5, 0])

