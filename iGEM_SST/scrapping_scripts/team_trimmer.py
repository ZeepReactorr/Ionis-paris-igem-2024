import os 
import pandas as pd 
from selenium import webdriver 
from selenium.webdriver import Chrome 
from selenium.webdriver.common.by import By 

os.chdir(BASE_DIR = '\\'.join(os.path.dirname(os.path.abspath(__file__)).split('\\')[:-1]))

all_team = open('finish.txt', 'r')
all_team_data = open('last100.txt', 'w', encoding='utf-8')
for i in all_team.readlines():
    url = i.split('\t')[0]
    village = i.split('\t')[1].strip('\n')

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
    try :
        content = driver.find_element(By.CSS_SELECTOR, "p[class*='ui-pt-2 ui-text-sm'")
        content_text = str(content.text)
    except:
        content_text = ''
            
    all_team_data.write(url + '\t' + village + '\t' + content_text + '\n')
    
all_team_data.close()
