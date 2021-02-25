from mesa_geo import GeoAgent, AgentCreator, GeoSpace
from mesa.time import BaseScheduler
from mesa import datacollection
from mesa import Model
from scipy import stats
import math
import os
import sys
import warnings


#Shapely Imports
from shapely.geometry import Polygon, Point, LineString
import shapely


#Data Analysis
import geopandas as gpd
import pandas as pd
import numpy as np
import random
import csv

# Configuration Data and Files
import configparser

#Plot
import matplotlib.pyplot as plt


# Prefix for config data
#os.chdir(os.path.dirname(sys.path[0]))
config_file_path_prefix = './config/'


# parser viz config data
viz_ini_file = 'vizparams.ini'

parser_viz = configparser.ConfigParser()
parser_viz.read(config_file_path_prefix + viz_ini_file)

default_section = parser_viz['DEFAULT_PARAMS']


# parser disease config data

disease_params_ini = 'diseaseparams.ini'
parser_dis = configparser.ConfigParser()
parser_dis.read(config_file_path_prefix + disease_params_ini)
incubation = parser_dis['INCUBATION']


# NPI config data


npi_params_ini = 'NPI.ini'
parser_npi = configparser.ConfigParser()
parser_npi.read(config_file_path_prefix + npi_params_ini)






def generate_random(polygon):
    minx, miny, maxx, maxy = polygon.bounds
    while True:
        pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
        if polygon.contains(pnt):
            return pnt






def find_room_type(room_agents, room_type):
    """
    room_agents: a list containing room agents
    room_type: a valid string of room type: [None, 'restroom_grade_boys', 'lunch_room', 'classroom_grade',
       'restroom_all', 'restroom_grade_girls', 'restroom_KG',
       'classroom_KG', 'community_room', 'library',
       'restroom_special_education', 'restroom_faculty',
       'classroom_special_education', 'health_room', 'faculty_lounge',
       'classroom_preschool', 'restroom_preschool']
    """
    out = []
    for agent in room_agents:
        if agent.room_type == room_type:
            out.append(agent)
    return out





def load_map(file_path):
    '''
    This is specific to the current school layout at this point, should be modified later in the future
    assume the input school map file is sufficient
    '''
    school_geometry = gpd.read_file(file_path)

    # move second floor to map bottom
    p = Polygon([(900,-800), (900,-1100), (1650,-1100), ( 1650,-800)])
    sf = school_geometry.geometry.intersects(p)
    new_sf = school_geometry[sf].translate(yoff=-800)
    school_geometry.loc[sf,["geometry"]] = new_sf
    school_gdf = school_geometry.rename(columns={"object": "room_type"})


    # generate recess area
    grd_recess = Polygon([(750, -1000), (750, -710), (1615, -710), (1615, -1480), (1450, -1480), (1450, -1200),
                                                               (900, -1200), (900, -1000)])
    pkg_recess = Polygon([(430, -1400), (430, -1150), (660, -1150), (660, -1400)])



    school_gdf = school_gdf[['Id', 'room_type', 'geometry']]
    school_gdf.loc[len(school_gdf)] = [90000, 'recess_yard', grd_recess]
    school_gdf.loc[len(school_gdf)] = [90001, 'recess_yard', pkg_recess]
    return school_gdf




class Human(GeoAgent):

    # plot config

    marker = default_section['marker']
    colordict = {"healthy": default_section['healthy'], 'exposed': default_section['exposed'], 'infectious': default_section['infectious']}
    edgedict = {"healthy": default_section['healthy_edge'], 'exposed': default_section['exposed_edge'], 'infectious': default_section['infectious_edge']}
    sizedict = {"healthy": default_section['healthy_size'], 'exposed': default_section['exposed_size'], 'infectious': default_section['infectious_size']}




#     # dummy config for data collection
#     viral_load = None

    shape, loc, scale = (float(incubation['shape']), float(incubation['loc']), float(incubation['scale']))
    x = np.linspace(-10, 8, 19)
    infective_df = pd.DataFrame({'x': list(x), 'gamma': list(stats.gamma.pdf(x, a=shape, loc=loc, scale=scale))})
    
    ## TODO: Once Humans are in a separate .py file
    # i.e. each agent type has different file
    ## 

#     # infectious curve
#     range_data = list(range(int(incubation['lower_bound']), int(incubation['upper_bound']) + 1))
#     infective_df = pd.DataFrame(
#         {'x': range_data,
#          'gamma': list(stats.gamma.pdf(range_data, a=shape, loc=loc, scale=scale))
#         }
#     )


    def __init__(self, unique_id, model, shape, room, health_status = 'healthy'):
        super().__init__(unique_id, model, shape)



        # disease config

        self.health_status = health_status
        prevalence = float(parser_dis['ASYMPTOMATIC_PREVALENCE']['prevalence'])
        self.asymptomatic = np.random.choice([True, False], p = [prevalence, 1-prevalence])
        self.symptoms = False

        
        self.infective = False

        # symptom onset countdown config
        
        
        shape, loc, scale =  (0.6432659248014824, -0.07787673726582335, 4.2489459496009125)
        x = np.linspace(0, 17, 1000)
        
        # initialize base transmission rate for infected students
        if self.health_status != 'healthy':
            temp = int(np.round(stats.lognorm.rvs(shape, loc, scale, size=1)[0], 0))
            if temp > 17:
                temp = 17
            if temp < -10:
                temp = 10
            
            self.time_to_symptom = int(np.round(stats.lognorm.rvs(shape, loc, scale, size=1)[0], 0))
        # create infectiveness reference dataframe
        shape, loc, scale = (20.16693271833812, -12.132674385322815, 0.6322296057082886)


        # temp Mask Passage: 1 = no masks, .1 = cloth, .05 = N95
        self.mask_passage_prob = .1
        
        # From 10000 lognorm values in R
        countdown = parser_dis['COUNTDOWN']
        shape, loc, scale =  (float(countdown['shape']), float(countdown['loc']), float(countdown['scale']))

        lognormal_dist = stats.lognorm.rvs(shape, loc, scale, size=1)

        num_days = min(np.round(lognormal_dist, 0)[0], int(countdown['upper_bound'])) # failsafe to avoid index overflow
        self.symptom_countdown = int(num_days)
        
        # room setup
        self.room = room
        self.x = self.shape.x
        self.y = self.shape.y






    def update_shape(self, new_shape):
        self.shape = new_shape
        self.x = self.shape.x
        self.y = self.shape.y


    def __update(self):
        '''
        # treat every 'exposed' student as infectious, then calculate likelihood of infection using time_until_symptoms
        
        stepwise infection calculator for each unique individual in the simulation
        
        1. baseline transmission generated from temporal analysis papers ** cite
        
        2. distance multiplier calculated from Chu paper ** cite
        
        3. mask passage probability justified by MIT / Cambridge (who also use .1 for cloth and .05 for N95)
        
        4. time steps are 5 minutes: multiply by 5/60 to get transmission rate/5 minute steps rather than transmission rate/hour
        
        final infectivity calculation: 1 * 2 * 3 * 4   
        
        
        '''
        # infection above max range is considered as aerosal transmission
        if self.mask and not (self.model.activity[self.room.schedule_id] == 'lunch'):
            neighbors = self.model.grid.get_neighbors_within_distance(self, int(parser_npi['MASKS']['infection_distance']))
        else:
            neighbors = self.model.grid.get_neighbors_within_distance(self, int(parser_npi['NO_NPI']['infection_distance']))
        



        if self.health_status == 'exposed' and self.infective: # if exposed

            # normalize symptom countdown value to infectious distribution value
            # 0 being most infectious
            # either -10 or 8 is proven to be too small of a chance to infect others, thus covering asympotmatic case


            countdown_norm = min(int(incubation['upper_bound']), max(int(incubation['lower_bound']), 0 - self.symptom_countdown))
            temp_prob = self.infective_df[self.infective_df['x'] == countdown_norm]['gamma'].iloc[0]


            for neighbor in neighbors: # could we do my inf/uninf dicts?

                # Check class is Human
                if issubclass(type(neighbor), Human):
                    
                    if (neighbor.unique_id != self.unique_id) and (neighbor.health_status !-: # and is not infected infecting another infected
                        # Call Droplet transmission function
                        temp_prob = droplet_infect(self, neighbor)
                                                                   
                        infective_prob = np.random.choice ([True, False], p = [temp_prob, 1-temp_prob])
                        if infective_prob and self.__check_same_room(neighbor):
                            neighbor.health_status = 'exposed'
                            
                            
                   
    def droplet_infect(infected, uninfected):
        '''
        baseline transmission rate
        '''
        
        distance = infected.shape.distance(uninfected.shape)
        time = infected.time_to_symptom
        transmission_baseline = infective_df[infective_df.x == -1 * time]['gamma']
                                                                   
     
        # Use Chu distance calculation ## see docs
        chu_distance_multiplier = 1/2.02
        distance_multiplier = distance * chu_distance_multiplier                                                           

        
        # approximate student time spent breathing vs talking vs loudly talking
        breathing_type_multiplier = np.random.choice([.1, .5, 1], p=[.2, .05, .75])
        # whisper, loud, heavy

        # temp Mask Passage: 1 = no masks, .1 = cloth, .05 = N95
        mask_passage_prob = .1
        mask_multiplier = mask_passage_prob # equivalent to cloth masks

        # convert transmission rate / hour into transmission rate / step
        hour_to_fivemin_step = 5/60

        return transmission_baseline * distance_multiplier * breathing_type_multiplier * mask_multiplier * hour_to_fivemin_step




    def __check_same_room(self, other_agent):
        '''
        check if current agent and other agent is in the same room

        the purpose of this function is to make sure to eliminate edge cases that one agent near the wall of its room
        infects another agent in the neighboring room

        this is at this iteration of code only implemented for class purpose, as unique id check is way more efficient

        later implementation should add attribute to human agent for current room

            other_agent: other agent to check
            returns: boolean value for if the two agents are in the same room
        '''
        same_room = True
        if self.model.activity[self.room.schedule_id] == 'class':
            same_room = (self.room.unique_id == other_agent.room.unique_id)
        return same_room


    def __move(self, move_spread = 12, location = None):
        '''
        Checks the current location and the surrounding environment to generate a feasbile range of destination (the area
        of a circle) for agent to move to.
        The radius of the circle is determined by agent's move_factor.
        Assigns new point to override current point.
        '''

        if not location:
            location = self.room
        move_spread = location.shape.intersection(self.shape.buffer(move_spread))
        minx, miny, maxx, maxy = move_spread.bounds
        while True:
            pnt = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
            # check if point lies in true area of polygon
            if move_spread.contains(pnt):
                self.update_shape(pnt)
                break


    def plot(self):
        plt.plot(
            self.shape.x, self.shape.y,
            marker=self.marker,
            mec = self.edgedict[self.health_status],
            color = self.colordict[self.health_status],
            markersize = self.sizedict[self.health_status]
                )











class Student(Human):
    def __init__(self, unique_id, model, shape, room, health_status = 'healthy', mask_on=False):
        super().__init__(unique_id, model, shape, room, health_status)

        viz_ini_file = 'vizparams.ini'
        parser = configparser.ConfigParser()
        parser.read(config_file_path_prefix + viz_ini_file)
        student_viz_params = parser['STUDENT']

        self.grade = self.room.room_type.replace('classroom_', '')
        self.mask = mask_on
        self.seat = Point(self.shape.x, self.shape.y)
        self.marker = student_viz_params['marker']




        self.out_of_place = False
        self.prev_activity = None
        self.lunch_count = 0











    def step(self):
        self._Human__update()

        activity = self.model.activity[self.room.schedule_id]
        # case 1: student in class
        if activity == 'class':
            if self.prev_activity != activity:
                self.prev_activity = activity
                self.update_shape(self.seat)

            if self.room.prob_move:
                self.out_of_place = True
                self._Human__move()
            else:
                if self.out_of_place:
                    self.update_shape(self.seat)
                    self.out_of_place = False



        # case 2: student in recess
        elif activity == 'recess':
            ### TODO: recess yard assignment is hard coded ###

            location = self.model.recess_yards[0]
            if self.grade != 'grade':
                location = self.model.recess_yards[1]

            if self.prev_activity != activity:
                self.update_shape(generate_random(location.shape))
                self.prev_activity = activity

            self._Human__move(move_spread=15, location = location)



        # case 3: student having lunch
        elif activity == 'lunch':
            #in class lunch case
            if self.model.inclass_lunch or self.grade != 'grade':
                if self.prev_activity != activity:
                    self.update_shape(self.seat)
                    self.prev_activity = activity
                    self.out_of_place = True
                    self._Human__move()
                else:
                    if self.out_of_place:
                        self.update_shape(self.seat)
                        self.out_of_place = False


            #in cafeteria lunch case
            else:
                if self.prev_activity != activity:
                    self.update_shape(generate_random(self.model.lunchroom.shape))
                    self.prev_activity = activity

                # enter lunch cafeteria, move free for 2 iteration
                if self.lunch_count < 2:
                    self._Human__move(move_spread=40, location = self.model.lunchroom)

                # finds seat, settle in seat until finish lunch
                elif self.lunch_count == 2:
                    self.update_shape(self.model.lunchroom.seats[0])
                    # remove seat from seat list if the student occupies it
                    self.model.lunchroom.seats = self.model.lunchroom.seats[1:]

                # release seat back to seat list
                elif self.lunch_count == 7:
                    self.model.lunchroom.seats.append(self.shape)
                    self.lunch_count = -1


                self.lunch_count += 1












class Teacher(Human):
    def __init__(self, unique_id, model, shape, room, health_status = 'healthy',mask_on=True):
        super().__init__(unique_id, model, shape, room, health_status)
        self.mask = mask_on
        self.classroom = self.room # TODO: for future development enabling Teachers to move to other room during non-class time

        viz_ini_file = 'vizparams.ini'
        parser = configparser.ConfigParser()
        parser.read(config_file_path_prefix + viz_ini_file)
        teacher_viz_params = parser['TEACHER']

        self.marker = teacher_viz_params['marker']
        self.edgedict = {"healthy": teacher_viz_params['healthy_edge'], 'exposed': teacher_viz_params['exposed_edge'], 'infectious': teacher_viz_params['infectious_edge']}
        self.sizedict = {"healthy": teacher_viz_params['healthy_size'], 'exposed': teacher_viz_params['exposed_size'], 'infectious': teacher_viz_params['infectious_size']}


    def step(self):
        self._Human__update()
        self._Human__move()







class Classroom(GeoAgent):


    # dummy config for data collection
    health_status = None
    symptoms = None
    x = None
    y = None


    def __init__(self, unique_id, model, shape, room_type):
        super().__init__(unique_id, model, shape)
        #self.occupants = occupants #List of all occupants in a classroom
        #self.barrier = barrier_type
        self.room_type = room_type
        self.viral_load = 0
        self.prob_move = False
        self.schedule_id = None

    def step(self):
        # roll weighted probability for current step having in-class activity
        self.prob_move = np.random.choice([True, False], p = [0.2, 0.8])


        # UPDATE 10/16: move aerosal update entirely to room object, better for future room volume inclusion
        # TODO: this aerosal version is very immature!

        occupants = [a for a in list(self.model.grid.get_intersecting_agents(self)) if issubclass(type(a), Human)]
        # this upates room viral_load to an arbitrary that's scaled with num of covid patients and area of room
        # should change to volume of room later
        self.viral_load += len([a for a in occupants if a.health_status != "healthy"])/self.shape.area

        # Don't run anything at this point
        # intended for aerosal transmission
        '''
        for agent in occupants:
            ### TODO: check literature for aerosal infection ###
            if issubclass(type(agent), Human):
                prob_exposed = self.viral_load/100

                if agent.mask:
                    prob_exposed -= 0.1 # particularly here mask effect on room aerosal transmission
                    if prob_exposed < 0:
                        prob_exposed = 0

                if np.random.choice([True, False], p = [prob_exposed, 1-prob_exposed]):
                    agent.health_status = 'exposed'
        '''



        # ventilation natrually
        ### TODO: check literature for aerosal infection ###
        self.viral_load = max(self.viral_load-0.002, 0)

        if self.schedule_id is not None:
            if "classroom" in self.room_type:
                if self.model.activity[self.schedule_id] != 'class':
                    self.ventilate()
        elif (self.room_type == 'lunch_room') and (not ('lunch' in self.model.activity.to_numpy())):
            self.ventilate()


    def ventilate(self):
        ### TODO: check literature for aerosal infection ###
        self.viral_load = max(self.viral_load-0.05, 0)


    def generate_seats(self, N, width):
        self.seats = []
        center = self.shape.centroid
        md = math.ceil(N**(1/2))
        pnt = Point(center.x - width*md//2, center.y - width*md//2)
        for i in range(md):
            for j in range(md+1):
                self.seats.append(Point(pnt.x + i*width, pnt.y + j*width))


    def generate_seats_lunch(self, xwidth, ywidth):

        self.seats = []
        xmin, ymin, xmax, ymax = self.shape.bounds
        xcoords = xmin + xwidth
        ycoords = ymin + ywidth

        y_pointer = ycoords
        x_pointer = xcoords

        while (xcoords < xmax):

            while (ycoords < ymax):
                self.seats.append(Point(xcoords, ycoords))
                ycoords += ywidth

            xcoords += xwidth
            ycoords = y_pointer

        np.random.shuffle(self.seats)





class School(Model):


    def __init__(self, map_path, schedule_path, grade_N, KG_N, preschool_N, special_education_N,
                 faculty_N, seat_dist, init_patient=3, attend_rate=1, mask_prob=0.516, inclass_lunch=False, username="jleiucsd"):
        # zipcode etc, for access of more realistic population from KG perhaps



        # model param init
        self.__mask_prob = mask_prob
        self.inclass_lunch = inclass_lunch
        self.seat_dist = math.ceil(seat_dist/(attend_rate**(1/2)))
        self.idle_teachers = [] # teachers to be assigned without a classroom
        self.init_patient = init_patient





        # mesa model init
        self.running = True
        self.grid = GeoSpace()
        self.schedule = BaseScheduler(self)



        #data collect init
        model_reporters = {"day": "day_count",
                           "cov_positive": "infected_count"}
        agent_reporters = {"unique_id": "unique_id",
                           "health_status": "health_status",
                           "symptoms": "symptoms",
                           "x": "x",
                           "y": "y",
                           "viral_load": "viral_load"}
        self.datacollector = datacollection.DataCollector(model_reporters=model_reporters, agent_reporters=agent_reporters)






        school_gdf = load_map(map_path)





        # room agent init
        self.room_agents = school_gdf.apply(lambda x: Classroom(
        unique_id = x["Id"],
        model = self,
        shape = x["geometry"],
        room_type = x["room_type"]),
                     axis=1
                    ).tolist()

        self.grid.add_agents(self.room_agents)



        # stats tracking init
        self.infected_count = 0
        self.step_count = 0
        self.day_count = 0
        self.num_exposed = 0



        # student activity init
        self.schoolday_schedule = pd.read_csv(schedule_path)
        self.activity = None


        # id tracking init
        self.__teacher_id = 0
        self.__student_id = 0
        self.__faculty_N = faculty_N
        self.schedule_ids = self.schoolday_schedule.columns



        self.recess_yards = find_room_type(self.room_agents, 'recess_yard')


        def init_agents(room_type, N, partition=False):
            '''
            batch initialize human agents into input room type rooms with equal partition size

            room_type: a valid string of room type: [None, 'restroom_grade_boys', 'lunch_room', 'classroom_grade',
               'restroom_all', 'restroom_grade_girls', 'restroom_KG',
               'classroom_KG', 'community_room', 'library',
               'restroom_special_education', 'restroom_faculty',
               'classroom_special_education', 'health_room', 'faculty_lounge',
               'classroom_preschool', 'restroom_preschool']
            '''


            rooms = find_room_type(self.room_agents, room_type)


            # if student group should be seperated to different day schedules
            # assigning schedule_id to equally partitioned rooms
            # currently only grade 1-5 "grade" students need to be partitioned,
            partition_size = len(rooms)
            if partition:
                partition_size = math.ceil(partition_size/len(self.schedule_ids))

            class_size = N//len(rooms)
            remaining_size = N%len(rooms)

            for i, classroom in zip(range(len(rooms)), rooms):

                classroom.generate_seats(class_size, self.seat_dist)
                classroom.schedule_id = self.schedule_ids[i//partition_size]

                for idx in range(class_size):
                    pnt = classroom.seats[idx]
                    mask_on = np.random.choice([True, False], p=[mask_prob, 1-mask_prob])
                    agent_point = Student(model=self, shape=pnt, unique_id="S"+str(self.__student_id), room=classroom, mask_on=mask_on)

                    self.grid.add_agents(agent_point)
                    self.schedule.add(agent_point)
                    self.__student_id += 1

                # spread remaining student into all classrooms
                if remaining_size > 0:
                    pnt = classroom.seats[class_size]
                    mask_on = np.random.choice([True, False], p=[mask_prob, 1-mask_prob])
                    agent_point = Student(model=self, shape=pnt, unique_id="S"+str(self.__student_id), room=classroom, mask_on=mask_on)

                    self.grid.add_agents(agent_point)
                    self.schedule.add(agent_point)
                    self.__student_id += 1
                    remaining_size -= 1


                #add teacher to class
                pnt = generate_random(classroom.shape)
                agent_point = Teacher(model=self, shape=pnt, unique_id="T"+str(self.__teacher_id), room=classroom)
                self.grid.add_agents(agent_point)
                self.schedule.add(agent_point)
                self.idle_teachers.append(agent_point)
                self.__teacher_id += 1
                self.__faculty_N -= 1



        # initialize all students and teachers in classrooms
        init_agents("classroom_grade", int(grade_N*attend_rate), partition=True)
        # keep track of student types
        #self.grade_students = [a for a in list(self.schedule.agents) if isinstance(a, Student)]
        init_agents("classroom_KG", int(KG_N*attend_rate))
        init_agents("classroom_preschool", int(preschool_N*attend_rate))
        #self.pkg_students = [a for a in list(set(self.schedule.agents).difference(self.grade_students)) if isinstance(a, Student)]
        init_agents("classroom_special_education", int(special_education_N*attend_rate))



        # dump remaining teacher to faculty lounge
        for f_lounge in find_room_type(self.room_agents, "faculty_lounge"):
            f_lounge.schedule_id = self.schedule_ids[0]

            for i in range(self.__faculty_N):

                pnt = generate_random(f_lounge.shape)
                agent_point = Teacher(model=self, shape=pnt, unique_id="T" + str(self.__teacher_id), room=f_lounge)
                self.grid.add_agents(agent_point)
                self.schedule.add(agent_point)
                self.__teacher_id += 1

        #self.people = list(self.schedule.agents)

        # add rooms to scheduler at last
        for room in self.room_agents:
            self.schedule.add(room)




        self.lunchroom = find_room_type(self.room_agents, 'lunch_room')[0]
        self.lunchroom.generate_seats_lunch(3, 12)






    def small_step(self):
        self.schedule.step()
        self.grid._recreate_rtree()



    def add_N_patient(self, N):
        patients = random.sample([a for a in self.schedule.agents if isinstance(a, Student)], N)
        for p in patients:
            p.health_status = "exposed"
            p.asymptomatic = True
            p.infective = True


    def show(self):
        '''
        plot current step visualization
        deprecated since end of model visualization update
        '''

        # UPDATE 10/16: add deprecation warning
        message  = "this function is no longer used for performance issues, check output_image.py for end of model visualization"
        warnings.warn(message, DeprecationWarning)


        school_geometry = gpd.GeoSeries([a.shape for a in self.room_agents])
        school_map = gpd.GeoDataFrame({"viral_load" : [min(a.viral_load, 5) for a in self.room_agents]})
        school_map.geometry = school_geometry
        basemap = school_map.plot(column = "viral_load", cmap="Reds", alpha = 0.5, vmin = 0, vmax = 5)
        school_map.boundary.plot(ax = basemap, color='k', linewidth=0.2)

        list(map(lambda a: a.plot(), [a for a in self.schedule.agents if issubclass(type(a), Human)]))

        hour = 9 + self.step_count*5//60 # assume plot start at 9am
        minute = self.step_count*5%60
        plt.title("Iteration: Day {}, ".format(self.day_count + 1) + "%d:%02d" % (hour, minute), fontsize=30)



    def __update_day(self):
        '''
        update incubation time, reset viral_load, remove symptomatic agents, etc for end of day
        '''


        for a in self.schedule.agents[:]:
            if issubclass(type(a), Human):


                if a.symptoms:
                    # remove agent if symptom onset
                    if isinstance(a, Teacher):
                        # assign a new teacher to position
                        new_teacher = self.idle_teachers.pop()
                        new_teacher.shape = a.shape
                        new_teacher.room = a.room
                        new_teacher.classroom = a.classroom
                    self.schedule.remove(a)
                    self.grid.remove_agent(a)

                # UPDATE 10/16: infectious made obsolete, end of day update rework
                elif a.health_status == "exposed":
                    # UPDATE 10/17: update infective delay if agent is not infective by end of day
                    a.infective = True
                    a.symptom_countdown -= 1
                    # calculate when symptoms begin to show using 0-15 density
                    if a.symptom_countdown <= 0:
                        if a.symptom_countdown == 0:
                            self.infected_count += 1
                        # update model stat for total infected
                        # negative countdown means this agent is asymptomatic

                        if not a.asymptomatic:
                            # this is a really small chance, however possible
                            # set symtoms to true
                            # next day this agent will be removed from the model
                            a.symptoms = True

            else:
                # reset viral_load of room agents
                a.viral_load = 0

    def step(self):
        '''
        simulate a day with school day schedule
        '''
        if not self.schedule.steps:
            self.add_N_patient(self.init_patient)



        for i, row in self.schoolday_schedule.iterrows():
            self.activity = row
            self.datacollector.collect(self)
            self.schedule.step()
            self.grid._recreate_rtree()
            self.step_count += 1


        self.__update_day()
        self.grid._recreate_rtree()
        self.day_count += 1
        self.step_count = 0
