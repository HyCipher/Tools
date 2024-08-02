import csv

# 读取name文件中的数据，并存储在一个集合中
def read_name_file(filepath):
    with open(filepath, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        names = set(row[0].strip() for row in reader if len(row) > 0)
        # print(names)
    return names


# 读取CSV文件的第四列数据，并检查是否存在于name集合中
def check_names_in_csv(csv_filepath, names_filepath):
    names = read_name_file(names_filepath)

    with open(csv_filepath, mode='r', newline='') as csvfile:
        reader = csv.reader(csvfile,delimiter='\t')
        # 遍历每行数据
        for row in reader:
            # print(row)
            if len(row) >= 3:
                name_to_check = row[3].replace('ID=gene-', '')  # 第四列索引为3
                # print(name_to_check)
                if name_to_check not in names:
                    print(f"Not found: {name_to_check}")


def launch(variable):
    csv_filepath = f'/BiO/Live/rooter/Downloads/ortholog/flanking/clean/{variable}_flanking-cleansing.csv'
    names_filepath = f'/BiO/Live/rooter/Downloads/ortholog/longest_bank/longest/{variable}.longest.txt'
    check_names_in_csv(csv_filepath, names_filepath)
    
    
if __name__ == '__main__':
    list_file = '/BiO/Live/rooter/Downloads/ortholog/list.csv'
    with open(list_file, mode='r', newline='') as file:
        reader = csv.reader(file)
   
        # 遍历每行数据
        for row in reader:
            for variable in row:
                launch(variable)
