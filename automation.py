import os


def create_folder(path):
    try:
        os.system(f'mkdir {path}')
    except:
        print('Folder already exists')


if __name__ == '__main__':
    create_folder('C:\\Users\\User\\Desktop\\test')