# 기존 소요시간과 새로운 소요시간을 입력받아 개선율(%) 계산

old_time = [1.68, 2.48, 4.82, 5.56, 9.36, 2.07, 3.00, 4.16, 6.44, 8.34, 10.54, 13.28]
new_time = [1.52, 1.84, 2.28, 3.87, 4.48, 1.53, 2.31, 2.76, 3.96, 4.53, 5.29, 8.05, 9.17, 11.78]

for old, new in zip(old_time, new_time):
    improvement_rate = ((old - new) / old) * 100
    print(f"개선율: {improvement_rate:.2f}%")
