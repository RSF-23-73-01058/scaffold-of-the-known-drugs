from bs4 import BeautifulSoup


def main():
    # insert path to the image here
    path =
    svg_soup = BeautifulSoup(open(path), 'lxml')
    groups = svg_soup.find_all('g')
    groups_list = {int(group['id'][6:]) for group in groups if group['id'][6:] != ''}
    print(*[i for i in range(1, 522) if i not in groups_list], sep=' ')
    print(len(groups_list))


if __name__ == '__main__':
    main()