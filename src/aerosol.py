from math import exp

def transmission_aerosol_percent(length, width, height, num_occupants, num_infectious, num_masks, \
                                 airchanges=3, filters = 0, UV=0, event_duration=5, background_CO2=415, \
                                 temperature=30, humidity=65, pressure=0.95):
    '''
    This function derives its attributes from the Aerosol Transmission for COVID Spreadsheet.
    
    Name: estimate_probability
    
    Args:
        length: length of room in feet.
        width: width of room in feet.
        height: height of room in feet.
        num_occupants: number of occupants in room
        num_infectious: number of infectious individuals in room
        num_masks: number of people wearing masks in room
        airchanges: Number of airchanges per hour
        filters: If there are extra protectional measures such as HEPA filters etc.
        UV: UV Index: Defaulted to 0
        event_duration: Duration of event (5 minutes)
        background_CO2: Outdoor background CO2 concentration in ppm.
        temperature: Temperature in degrees celcius.
        humidity: Percentage humidity 
        pressure: Atmospheric pressure
        
    '''
    
    #Dimensions of room (converting to meters)
    length *= 0.305
    width *= 0.305
    height *= 0.305
    
    event_duration_hrs = event_duration / 60
    
    #Healthy people
    num_susceptible = num_occupants - num_infectious
    
    #Volume of room (in m^3)
    volume = length * width * height
    
    #Derived from sheet as provided by Dr. Paul Dabisch, Dept. of Homeland Security, USA. Units in h-1 
    decayrate = (7.56923714795655+1.41125518824508*(temperature-20.54)/10.66 + 0.02175703466389*(humidity-45.235)/28.665+7.55272292970083*((UV*0.185)-50) / 50 + (temperature-20.54)/10.66*(UV*0.185-50)/50*1.3973422174602)*60
    
    #Hardcoded in sheet
    surface_deposition = 0.3
    
    
    #Ventilation
    first_order_loss = decayrate + surface_deposition + airchanges + filters
    ventilation_rate_per_person = volume * (airchanges + filters) * 1000/3600
    ventilation_rate_per_person = ventilation_rate_per_person / num_occupants #In order to get units as L/s/person
    
    
    #Area and Density
    area = length * width
    area_per_person = area / (0.305**2)
    area_per_person = area_per_person / num_occupants #sq.ft / person
    
    density_room = num_occupants / area #persons / m^2
    density_room_volume = volume / num_occupants # m^3 / person
    
    
    #Activity Parameters
    breathing_rate = 0.0086*60
    
    #People determined Virus Parameters
    quanta_exalation_rate = 25 ######Uncertain: Acitvity and Person dependent#####
    exhalation_mask_efficiency = 0.5 
    inhalation_mask_efficiency = 0.3 
    mask_fraction = (num_masks / num_occupants) 
    
    infection_probability = 0.2 ####AREA/TIME SPECIFIC#######
    
    #####Results#######
    net_emission_rate = quanta_exalation_rate * (1 - exhalation_mask_efficiency*mask_fraction) * num_infectious
    avg_quanta_conc = net_emission_rate/first_order_loss/volume*(1-(1/first_order_loss/event_duration_hrs)*(1-exp(-first_order_loss*event_duration_hrs)))
    quanta_per_person = avg_quanta_conc * breathing_rate * event_duration_hrs * (1- inhalation_mask_efficiency * mask_fraction)
    prob_infection = 1 - exp(-quanta_per_person)
    
    
    return prob_infection