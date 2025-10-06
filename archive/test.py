from socket import timeout
import requests
from bs4 import BeautifulSoup

url = "https://casdc.china-vo.org/mirror/ALMAGAL/2019.1.00195.L/QA2"

res0 = requests.get(url, timeout=30)
res0.raise_for_status()
soup0 = BeautifulSoup(res0.content, 'html.parser')
links0 = soup0.find_all('a', href=True)

for link0 in links0:
    href0 = link0['href']
    if href0.endswith('/') and 'science_goal' in href0:
        sci_url = url.rstrip('/') + '/' + href0
        # print(href0)
        
        res1 = requests.get(sci_url, timeout=30)
        soup1 = BeautifulSoup(res1.content, 'html.parser')
        links1 = soup1.find_all('a', href=True)

        for link1 in links1:
            href1 = link1['href']
            if href1.endswith('/') and 'group.' in href1:
                group_url = sci_url.rstrip('/') + '/' + href1
                # print(href1)

                res2 = requests.get(group_url, timeout=30)
                soup2 = BeautifulSoup(res2.content, 'html.parser')
                links2 = soup2.find_all('a', href=True)

                for link2 in links2:
                    href2 = link2['href']
                    if href2.endswith('/') and 'member.' in href2:
                        mem_url = group_url.rstrip('/') + '/' + href2
                        # print(href2)

                        res3 = requests.get(mem_url, timeout=30)
                        soup3 = BeautifulSoup(res3.content, 'html.parser')
                        links3 = soup3.find_all('a', href=True)
                        # print(res3)
                        # print(soup3)
                        # print(links3)


                        for link3 in links3:
                            href3 = link3['href']
                            # print(href3)
                            if href3 == 'product/':
                                pro_url = mem_url.rstrip('/') + '/' + href3
                                # print(pro_url)

                                # break

                                res4 = requests.get(pro_url, timeout=30)
                                soup4 = BeautifulSoup(res4.content, 'html.parser')
                                links4 = soup4.find_all('a', href=True)

                                for link4 in links4:
                                    href4 = link4['href']
                                    if href4.endswith('cont.I.pbcor.fits'):
                                        fits_url = pro_url.rstrip('/') + '/' + href4
                                        print(fits_url)

