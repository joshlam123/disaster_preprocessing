from bs4 import BeautifulSoup
from collections import defaultdict
from tqdm import tqdm
from requests import get

import pandas as pd


class trawler():
    def find_item_in_url(soup, tag, seq=-1):
        cnt = 0 
        for k,v in tag.items():
            if cnt == 0:
                if k == "a":
                    if len(v) > 0:
                        res = soup.find_all(k, {"class":v}, href=True)
                    else:
                        res = soup.find_all(k, href=True)

                    res = [i['href'] for i in res]
                else:
                    if len(v) > 0:
                        res = soup.find_all(k, {"class": v})
                    else:
                        res = soup.find_all(k)
            else:
                if k == "a":
                    if len(v) > 0:
                        res = soup.find_all(k, {"class":v}, href=True)
                    else:
                        res = soup.find_all(k, href=True)
                    res = [i['href'] for i in res]
                    
                else:
                    if len(v) > 0:
                        res = res.find_all(k, {"class": v})
                    else:
                        res = res.find_all(k)
        
        if seq == "all":
            return res[:]
        else:
            return res[seq]

    def get_url(link):
        page = get(link)
        soup = BeautifulSoup(page.text, 'html.parser')
        return soup

    def fetch_details(list_of_tags, ctry = "Indonesia", region = "jakarta"):
        initialurl=f"http://adinet.ahacentre.org/reports/fetch_reports?page=1&c%5B%5D=6&cff%5B0%5D%5B%5D=1&cff%5B0%5D%5B%5D={region}&cff%5B1%5D%5B%5D=24&cff%5B1%5D%5B%5D={ctry}"
        # last_pg_xpath = "/html/body/div[4]/table/tbody/tr/td[2]/ul/li[9]/span/a"
        
        soup = get_url(initialurl)

        corresp_item = find_item_in_url(soup, list_of_tags)
        
        iter_len = corresp_item.get_text()

        all_items = list()

        for entry in range(1, int(iter_len)+1):
            print(entry)
            filterurl=f"http://adinet.ahacentre.org/reports/fetch_reports?page={entry}&c%5B%5D=6&cff%5B0%5D%5B%5D=1&cff%5B0%5D%5B%5D={region}&cff%5B1%5D%5B%5D=24&cff%5B1%5D%5B%5D={ctry}"
            soup = get_url(filterurl)

            list_of_filter_tags = {"div":"rb_report verified", "a":"r_title"}

            filtered_items = find_item_in_url(soup, list_of_filter_tags, "all")
        
            all_items.append(filtered_items)

        return all_items
            #rb_list-view


class handleText():
    def stripProperly(txt):
        return txt.replace('\n',' ').replace('\t',' ').replace('\r','')

    def handleException(result):
        if result:
            return result.text
        else:
            return "" 

    def stripDescription(pageSoup):
        souped = pageSoup.find('div', class_='report-description-text')
        if souped is None:
            return "", "", "", "", ""
        else:
            souped = souped.get_text().strip('')

            desc = stripProperly(souped[12:souped.find('IMPACT')].strip())
            impact = stripProperly(souped[souped.find("IMPACT")+13:souped.find("RESPONSE")].replace('-',' ')).strip()
            response = stripProperly(souped[souped.find("RESPONSE")+13:souped.find("Additional Data")].replace('-',' ')).strip()
            additional = stripProperly(souped[souped.find("Additional Data")+16:souped.find("CasualtiesDamages")].replace('-',' ')).strip()
            casualties = stripProperly(souped[souped.find("CasualtiesDamages")+len('CasualtiesDamages')+1:].replace('-',' ')).strip()
        return desc, impact, response, additional, casualties

    def replaceOperation(soup):
        return soup.find_all('div', class_='report-category-list')[0].get_text().replace('\n','').replace('\t','').replace('\r','').replace('\xa0','')

    def getDetails(pageSoup):
        desc, impact, response, additional, casualties = stripDescription(pageSoup)
        if desc != "":
            pageDict = {'title':pageSoup.find_all('h1', class_='report-title')[0].get_text(), 'date' : pageSoup.find_all('span', class_='r_date')[0].get_text(),
                        'location' : pageSoup.find_all('span', class_='r_location')[0].get_text(), 'category' : replaceOperation(pageSoup),
                        'description' : desc, "impact" : impact, 'response' : response, 'additional' : additional, 'casualties' : casualties}
            return pageDict
        return

def get_all_url():
    tag, inner_tag, text  = "ul", "li", "span"
    list_of_tags = {tag:"pager", inner_tag:"", text:""} # to be traversed sequentially

    all_items = fetch_details(list_of_tags, "Indonesia", "jakarta")

    all_urls = list()
    for item in all_items:
        all_urls.extend(item)

    return all_urls, all_items

def create_masterlist(all_items):
    masterlist = defaultdict()
    urls = list()
    for item in all_items:
        urls.extend(item[:])

    for idx, url in enumerate(urls):
        page = requests.get(url)
        soup = BeautifulSoup(page.text, 'html.parser')
        details = getDetails(soup)
        if details is not None:
            masterlist[idx] = details

    for k,v in masterlist.items():
        location = v['title'].split(',')[0]
        masterlist[k]['location'] = location

    return masterlist

def make_masterlist_data(masterlist, all_urls):
    df = pd.DataFrame(masterlist).T
    df['url'] = all_urls
    return df
    

def main(save_dir):
    all_urls, all_items = get_all_url()
    masterlist = create_masterlist(all_items)
    df = make_masterlist_data(masterlist, all_urls)

    # do some prechecks here.
    df.to_csv("ADINET_indo.csv")

    print("Successfully saved ")

if __name__ == '__main__':
    # get input country and all that here. 
    main()


