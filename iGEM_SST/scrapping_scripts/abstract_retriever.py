from selenium import webdriver 
from selenium.webdriver import Chrome 
from selenium.webdriver.common.by import By 
import os, sys

os.chdir(BASE_DIR = '\\'.join(os.path.dirname(os.path.abspath(__file__)).split('\\')[:-1]))


def igem_sniffer(teams_file, words_test):
    teams = open(teams_file, 'r')
    
    to_consider = open('to_consider.txt', 'w')
    for i in teams.readlines():
        url = i.strip('\n')
        print(url)
        # Define the Chrome webdriver options
        options = webdriver.ChromeOptions() 
        options.add_argument("--headless") # Set the Chrome webdriver to run in headless mode for scalability

        # By default, Selenium waits for all resources to download before taking actions.
        # However, we don't need it as the page is populated with dynamically generated JavaScript code.
        options.page_load_strategy = "none"

        # Pass the defined options objects to initialize the web driver 
        driver = Chrome(options=options) 
        # Set an implicit wait of 5 seconds to allow time for elements to appear before throwing an exception
        driver.implicitly_wait(5)
         
        driver.get(url) 
        content = driver.find_element(By.CSS_SELECTOR, "div[class*='ui-content-box__text ui-h-full'")
        content_text = content.text
        
        for i in words_test:
            if i.lower() in str(content_text).lower():
                to_consider.write(url + '\n')     

    teams.close()
    to_consider.close()    
    
if __name__ == "__main__":
    datafile, keywords_list = sys.argv[1], sys.argv[2]
    igem_sniffer('liste_project_bioremédiation.txt', ['algae'])
