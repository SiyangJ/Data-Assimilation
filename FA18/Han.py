import os

if __name__ == '__main__':
    a_rows = int(input().strip())
    a_columns = int(input().strip())
    a = []
    for _ in range(a_rows):
        a.append(list(map(int, input().rstrip().split())))
    print(a[0][0],type(a[0][0]))
    print(a)
    