from mesa_geo import GeoAgent, AgentCreator, GeoSpace
from mesa.time import BaseScheduler
from mesa import datacollection
from mesa import Model
from scipy import stats 
import math
import os


#Shapely Imports
from shapely.geometry import Polygon, Point, LineString
import shapely


#Data Analysis 
import geopandas as gpd
import pandas as pd
import numpy as np
import random
import csv


#Plot
import matplotlib.pyplot as plt
import gc






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
    marker = 'o'
    colordict = {"healthy": '#00d420', 'exposed': '#f0e10c', 'infectious': '#a82e05'}
    edgedict = {"healthy": '#00d420', 'exposed': '#cf9e19', 'infectious': '#a82e05'}
    sizedict = {"healthy": 5, 'exposed': 6, 'infectious': 8}
    
    # disease config
    asymptomatic = np.random.choice([True, False], p = [0.58/100, 1- (0.58/100)])
    symptoms = False
    
    
    #dummy config for data collection
    viral_load = None
    
    def __init__(self, unique_id, model, shape, room, health_status = 'healthy'):
        super().__init__(unique_id, model, shape)
        self.health_status = health_status
        self.incubation_time = stats.gamma(5.807, 0.948).rvs(1)[0]
        self.room = room
        self.x = self.shape.x
        self.y = self.shape.y


        
        

 
    def update_shape(self, new_shape):
        self.shape = new_shape
        self.x = self.shape.x
        self.y = self.shape.y
    
    def __update(self):
        
        
        if self.health_status == 'infectious':
            # update room aerosal viral concentration:
            # Checks if infectious agent's neighbors are within 3 
            neighbors = self.model.grid.get_neighbors_within_distance(self, 12)
            if self.mask and not (self.model.activity[self.room.schedule_id] == 'lunch'):
                neighbors = self.model.grid.get_neighbors_within_distance(self, 4)

            # Loop through infected agent's neighbors that are within 3
            for neighbor in neighbors:

                # Check class is Human                            
                if issubclass(type(neighbor), Human):
                    protection_prob = True
                    if neighbor.mask:
                        protection_prob = np.random.choice ([True, False], p = [0.857, 1-0.857])
                        if self.model.activity[self.room.schedule_id] == 'lunch':
                            protection_prob = True


                    same_room = True
                    
                    # wall hard blocking transmission should only be an edge case if agents are in classroom
                    if self.model.activity[self.room.schedule_id] == 'class':
                        same_room = (self.room.unique_id == neighbor.room.unique_id)
                    if same_room and (neighbor.health_status == 'healthy'):
                        prob_exposed = np.random.choice ([True, False], p = [0.5, 0.5])                                    
                        if prob_exposed and protection_prob:
                            neighbor.health_status = 'exposed'
 
            
    
            
        if self.health_status == "exposed":
            self.incubation_time -= 0.00347222

            # modified infection based on step function update in student class
            if self.incubation_time <= 2:

                
                # TODO: literature proof for chance of testing positive ###
                prob_infect = np.random.choice([True, False], p=[0.8, 0.2])
                if prob_infect:
                    self.health_status = "infectious"
                else:
                    self.health_status = "healthy"
                    self.incubation_time = stats.gamma(5.807, 0.948).rvs(1)[0]


        if self.health_status == "infectious" and not self.symptoms:
            self.incubation_time -= 0.00347222
            if self.incubation_time <= 0:
                self.model.infected_count += 1    
                self.symptoms = True
    
    
    
    
    def __move(self, move_spread = 12, location = None):
        '''
        Checks the current location and the surrounding environment to generate a feasbile range of destination (the area
        of a circle) for agent to move to.
        The radius of the circle is determined by agent's move_factor.
        Assigns new point to override current point.
        '''   
        # a scale for health status
        if not location:
            location = self.room
        move_spread = location.shape.intersection(self.shape.buffer(40))
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

        self.grade = self.room.room_type.replace('classroom_', '')
        self.mask = mask_on        
        self.seat = Point(self.shape.x, self.shape.y)
        self.marker = 'o'
        
        
        

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
                    
            # update room aerosal virus concentration 
            if self.health_status == 'infectious':
                if self.mask:
                    self.room.viral_load += (0.01*0.01)
                else:
                    self.room.viral_load += 0.01
                
                
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
                # update room aerosal virus concentration 
                if self.health_status == 'infectious':
                    self.room.viral_load += 0.01

                    
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

                    
                if self.health_status == 'infectious':
                    self.model.lunchroom.viral_load += 0.01
            
                self.lunch_count += 1
    
                        
                            
                            
                            




                
                
    
class Teacher(Human):
    def __init__(self, unique_id, model, shape, room, health_status = 'healthy',mask_on=True):
        super().__init__(unique_id, model, shape, room, health_status)
        self.mask = mask_on
        self.classroom = self.room # TODO: for future development enabling Teachers to move to other room during non-class time
        
        self.marker = "^"
        self.edgedict = {"healthy": '#009416', 'exposed': '#cf9e19', 'infectious': '#a82e05'}
        self.sizedict = {"healthy": 7, 'exposed': 8, 'infectious': 10}
    
    
    def step(self):
        self._Human__update()
        self._Human__move()
        # update room aerosal virus concentration 
        if self.health_status == 'infectious':
            if self.mask:
                self.room.viral_load += (0.01*0.01)
            else:
                self.room.viral_load += 0.01
                
        





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
        self.prob_move = np.random.choice([True, False], p = [0.2, 0.8])
        for agent in self.model.grid.get_intersecting_agents(self):
            ### TODO: check literature for aerosal infection ###
            if issubclass(type(agent), Human):
                prob_exposed = self.viral_load/100
                
                if agent.mask:
                    prob_exposed -= 0.1 # particularly here mask effect on room aerosal transmission
                    if prob_exposed < 0: 
                        prob_exposed = 0
                
                if np.random.choice([True, False], p = [prob_exposed, 1-prob_exposed]):
                    agent.health_status = 'exposed'
        
        # ventilation natrually
        ### TODO: check literature for aerosal infection ###
        self.viral_load = max(self.viral_load-0.02, 0)
        
        if self.schedule_id is not None:
            if "classroom" in self.room_type:
                if self.model.activity[self.schedule_id] != 'class':
                    self.ventilate()
        elif (self.room_type == 'lunch_room') and (not ('lunch' in self.model.activity.to_numpy())):
            self.ventilate()
            
        
    def ventilate(self):
        ### TODO: check literature for aerosal infection ###
        self.viral_load = max(self.viral_load-0.5, 0)
        
        
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
        
        # model param init
        self.__mask_prob = mask_prob
        self.inclass_lunch = inclass_lunch
        self.seat_dist = math.ceil(seat_dist/(attend_rate**(1/2)))
        self.idle_teachers = [] # teachers to be assigned without a classroom
        self.init_patient = init_patient

        
        
        # output path init
        '''
        self.path = "/oasis/scratch/comet/{}/temp_project/".format(username) +\
        "output" +\
        "_attendrate_" + str(attend_rate) +\
        "_maskprob_" + str(self.__mask_prob) +\
        "_inclasslunch_" + str(self.inclass_lunch)
        
        try:
            os.mkdir(self.path)
        except OSError:
            pass
        '''

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
            p.health_status = "infectious"
            p.asymptomatic = True
        self.infected_count += N
    
    
    def show(self):

        ### TODO: update viral load to based on literature!! ###
        school_geometry = gpd.GeoSeries([a.shape for a in self.room_agents])
        school_map = gpd.GeoDataFrame({"viral_load" : [min(a.viral_load, 5) for a in self.room_agents]})
        school_map.geometry = school_geometry
        basemap = school_map.plot(column = "viral_load", cmap="Reds", alpha = 0.5, vmin = 0, vmax = 5)
        school_map.boundary.plot(ax = basemap, color='k', linewidth=0.2)

        list(map(lambda a: a.plot(), [a for a in self.schedule.agents if issubclass(type(a), Human)]))
        
        hour = 9 + self.step_count*5//60 # assume plot start at 9am
        minute = self.step_count*5%60
        plt.title("Iteration: Day {}, ".format(self.day_count + 1) + "%d:%02d" % (hour, minute), fontsize=30)
        
        '''
        img_path = self.path + "/outputimages"
        try:
            os.mkdir(img_path)
        except OSError:
            pass
        fileout = img_path + "/" + str(self.day_count + 1) + '_' + ("%03d" % (self.step_count)) + ".png"
        
        plt.savefig(fileout)
        # Clear the current axes.
        plt.cla() 
        # Clear the current figure.
        plt.clf() 
        # Closes all the figure windows.
        plt.close('all')
        # force clear catch 
        gc.collect()
        '''
        
    
    def __update_day(self):
        '''
        update incubation time, reset viral_load, remove symptomatic agents, etc for end of day
        '''
        steps_in_day = 60*24/5
        
        for a in self.schedule.agents[:]:
            if issubclass(type(a), Human):
                if a.health_status == "exposed":
                    a.incubation_time -= 0.00347222*(steps_in_day - self.step_count)
                     
                elif a.health_status == "infectious" and a.symptoms:
                    if isinstance(a, Teacher):
                        new_teacher = self.idle_teachers.pop()
                        new_teacher.shape = a.shape
                        new_teacher.room = a.room
                        new_teacher.classroom = a.classroom
                    self.schedule.remove(a)
                    self.grid.remove_agent(a)

            else:
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
