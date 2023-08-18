# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 15:44:00 2022

@author: Administrator
"""
from selenium import webdriver
from selenium.webdriver.common.by import By
import os
import re
import time

os.chdir('C:\\Users\\Administrator\\Documents\\Hettich\\SLT01')

def rupee_search(pdb_file):
    driver = webdriver.Firefox()
    driver.get("https://ayoubresearch.com/rupee/")
    
    pdb_radio = driver.find_element(By.XPATH, '/html/body/section[1]/div[2]/form/div[2]/div/label[2]/input')
    pdb_radio.click()
    
    contained_in = driver.find_element(By.XPATH, '/html/body/section[1]/div[2]/form/div[6]/div/select/option[2]')
    contained_in.click()
    
    file_path = driver.find_element(By.XPATH, '//*[@id="pdbFile"]')
    file_path.send_keys(pdb_file)
    
    search = driver.find_element(By.XPATH, '//*[@id="search"]')
    search.click()
    
    export = driver.find_element(By.XPATH, '/html/body/section[2]/div/div[2]/button')
    data_table = driver.find_element(By.XPATH, '//*[@id="dataTable"]')
    while not data_table.is_displayed():
        time.sleep(10)
    driver.execute_script("return arguments[0].scrollIntoView();", export)
    time.sleep(1)
    export.click()
    time.sleep(1)
        
    driver.quit()

    locus = re.search(r'PP_\d{4}', pdb_file).group()
    os.rename('C:\\Users\\Administrator\\Downloads\\data.csv', out_dir + locus + '-results.txt')
    
    time.sleep(15)

pdb_dir = os.path.abspath('alphafold/dwnld_unann') + '\\'
pdbs = [pdb_dir + f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
out_dir = os.path.abspath('alphafold/rupee_out') + '\\'

finished_files = set([re.search(r'PP_\d{4}', f).group() for f in os.listdir(out_dir)])
pdbs = [p for p in pdbs if not re.search(r'PP_\d{4}', p).group() in finished_files]


for pdb in pdbs:
    rupee_search(pdb)

