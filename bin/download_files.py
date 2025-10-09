'''
Download all pdf/ppt/pptx links in a url
'''
import os
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup

url = "https://cpraveen.github.io/talks.html"

#If there is no such folder, the script will create one automatically
folder_location = r'/tmp/download'
if not os.path.exists(folder_location):os.mkdir(folder_location)

response = requests.get(url)
soup = BeautifulSoup(response.text, "html.parser")     
pdf = soup.select("a[href$='.pdf']")
ppt = soup.select("a[href$='.ppt']")
pptx = soup.select("a[href$='.pptx']")
all_links = [pdf, ppt, pptx]

for links in all_links:
    for link in links:
        # Name pdf files using last portion of each link
        filename = os.path.join(folder_location,link['href'].split('/')[-1])
        print(filename)
        with open(filename, 'wb') as f:
            f.write(requests.get(urljoin(url,link['href'])).content)
