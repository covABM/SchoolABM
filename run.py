import sys
sys.path.append('src')


# data Analysis 
import geopandas as gpd
import pandas as pd
import numpy as np


# plot
import matplotlib.pyplot as plt
import gc

# output
from output_image import write_output

# multiprocessing
import multiprocessing
import brun

# model
from school_model import School

# spreadsheets
# import spreadsheets

#config
import configparser
import warnings
warnings.filterwarnings("ignore")


school_params = './config/schoolparams.ini'
parser = configparser.ConfigParser()
parser.read(school_params)

map_path = parser['SHAPEFILE']['path']
schedule_path = parser['SCHEDULE']['path'] 
schedule_steps = int(parser['TIME']['schedule_steps']) 


# if you wish to use google sheet parameter input, 
# you'll have to follow instructions for getting credentials of the sheet api
# follow the first two steps in:
# https://developers.google.com/sheets/api/quickstart/python
# and save your credentials in the src folder as listed below
# SHEET_URL = 'https://docs.google.com/spreadsheets/d/1Quyyey5B_kdQK1_OU0OkIDGZE27tuUfjZ6hCsV99-sM'
# credentials = './src/credentials.json'


####Pulling in parameters from config files####
population = parser['SCHOOL_POPULATION']
grade_N = int(population['grade_N'])
KG_N = int(population['KG_N'])
preschool_N = int(population['preschool_N'])
special_education_N = int(population['special_education_N'])
faculty_N = int(population['faculty_N'])

mask_prob = float(population['mask_prob'])
init_patient = int(population['init_patient'])
attend_rate = int(population['population_attending'])

days = int(parser['TIME']['days'])
max_steps = days*schedule_steps
iterations = 1

# This in the process of becoming a user adjustable parameter, 
# Hence, the alterative file path. 

activity_params = './config/activities.ini'
parser_activity = configparser.ConfigParser()
parser_activity.read(activity_params)
seat_dist = int(parser_activity['CLASSROOM_REGULAR']['seat_dist'])

school = School(map_path, schedule_path, grade_N, KG_N, preschool_N, special_education_N, 
                 faculty_N, seat_dist, init_patient=init_patient, attend_rate=attend_rate, mask_prob=mask_prob, inclass_lunch=True, username="jleiucsd")


while school.running and school.schedule.steps < max_steps:
    school.step()

agent_df = school.datacollector.get_agent_vars_dataframe()
model_df = school.datacollector.get_model_vars_dataframe()
print(agent_df.head())
print(model_df.head())
print(f'Model Simulation Complete for {days} days.')
