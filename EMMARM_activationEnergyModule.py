"""Estimation of apparent activation energy from EMM-ARM tests

This script allows the user to...

This file can also be imported as a module and contains the following
functions:
    * ...
    * ...
"""


import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib
matplotlib.use('QtAgg')
plt.style.use('tableau-colorblind10')
cm = 1/2.54  # centimeters in inches for plots
R=8.31446261815324/1000 #Universal gas constant


class EMMARMspecimen():
    """
    A class that defines all a EMMARM test specimen

    Attributes:
    - ... 

    Methods:
    - ...
    """

    def __init__(self, specimenID, dataAcquisitionSystem=None, samplingFrequency=None, systemChannel = 1):
        """Initialization of the class.

        Parameters
        ----------
        - specimenID: str
            Defines the identification used to refer to this specimen
        - dataAcquisitionSystem: str, optional
            Defines the acquisition system used to monitor this specimen.
            Options are: 'uEMMARM','National','old_uEMMARM','RPi'.
            If None is set, processing of acceleration data will not be possible for this specimen.
        - samplingFrequency: float, optional
            Defines the sampling frequency configured in the data acquisition system.
            If None is set, processing of acceleration data will not be possible for this specimen.
        - systemChannel: int, optional
            Defines what is the channel to be read from the file. Starts at 1 and unitarily increases (1,2,3...)
            If only one channel is available, then it is 1
        """

        ## Attributes defined in class instantiation
        self.specimenID = specimenID
        self.dataAcquisitionSystem = dataAcquisitionSystem
        #Validation of sampling frequency data based on the informed data acquistion system:
        if self.dataAcquisitionSystem == "uEMMARM":
            #uEMMARM is based on a fixed duration of measurement session, automatically obtained from test files
            self.samplingFrequency = "toBeDerivedFromExperimentalData"
        else:
            self.samplingFrequency = samplingFrequency
        self.systemChannel = systemChannel

        ## Attributes to be defined by other methods
        #Here they are set to None so a call of a method that depends on them, before a call of a method that sets them, can verify if their value is None and raise an error if so.
        
        #Attributes related to signal processing and modal analysis
        self.nps=None
        self.filterConfiguration=None
        #TODO: I think the below two attributes will not be maintained, because it is simpler to make methods of them and call when we want in the code
        self.modalIdentificationMethod_SingleAnalysis = {}
        self.modalIdentificationMethod_BatchAnalysis = {}
        
        #Attributes related to elastic modulus estimation from identified natural frequency values
        self.ageAtBeginningOfTest=None
        self.massOfEmptyTube=None
        self.totalLengthOfTube=None
        self.internalTubeDiameter=None
        self.externalTubeDiameter=None
        self.freeCantileverLengthEmptyTube=None
        self.frequencyEmptyTube=None
        self.massAtTipEmptyTube=None
        self.freeCantileverLength=None
        self.massOfFilledTube=None
        self.massAtTip=None

        #Atributes that will store values that are obtained after processign test data
        #There is no self.testTimeSeries_frequency since testNaturalFrequencySeries is directly mapped to testElasticModulusSeries
        self.testTimeSeries_elasticModulus=None 
        self.testNaturalFrequencySeries=None 
        self.testElasticModulusSeries=None
        #Temperature may not be monitored by the same system that monitors acceleration -> frequency -> elastic modulus
        #Therefore, its measurements may be related to different time series
        self.testTimeSeries_temperature=None 
        self.testTemperatureSeries=None

        #Attributes related to fitted model
        self.fittedEquationParameters=None
        #self.fittedStandardError stores one standard deviation errors associated to each parameter in the fitted model
        #Never use R2 in non linear fitting as fitting quality parameter
        #see: https://blog.minitab.com/en/adventures-in-statistics-2/why-is-there-no-r-squared-for-nonlinear-regression
        self.fittedStandardError=None
        self.fittedTimeSeries=None
        self.typeOfTimeSeriesUsedInFitting=None #Can be: "experimental_elasticModulus", "experimental_temperature", "custom"
        self.fittedElasticModulus=None
        self.fittedAdvancementDegree=None
        self.fittedDerivativeAdvancementDegree=None

        #Attributes to be used in the context of the Speed Method for calculating apparent activation energy
        self.fittedTimeSeries_EaCalculation=None
        self.derivativeAdvancementDegree_EaCalculation=None
        self.temperatureSeries_ForEaCalculation=None

        #Attributes set in the context of the equivalent time
        self.equivalentAge=None

        #Attributes related to reference methods that provides E-modulus estimates
        self.classicCompression_Time=None 
        self.classicCompression_Modulus=None 

    def readProcessedTestData(self, filePath, typeOfData):
        """Read file with "frequency", "elasticModulus", "temperature" data 

        A method to allow reading a .csv file which already contains processed test data.
        The proper attributes of the object will be populated with the file data.
        It may be used to read any property that is derived from the processing of raw experimental data: natural frequencies, elastic modulus, or temperature time series
        The file should be formatted in the following way:
            - The first line is the header containing the labels of each column (this is not read by the method)
            - 2 columns with a ";" separating the variables
            - The 1st column is the time series, in seconds
            - The 2nd column is the processed value (natural frequency, elastic modulus, or temperature)
            - Natural frequency should be in Hertz
            - Elastic modulus should be in GPa
            - Temperature should be in Celsius
        
        Parameters
        ----------
        - filePath: str
            Full path to the file
        - typeOfData: str
            The type of data to be read
            Options: "frequency", "elasticModulus", "temperature", "classicCompression"
        """

        if typeOfData not in [ "frequency", "elasticModulus", "temperature", "classicCompression"]:
            print("Type of data not compatible with this method. Specify either a 'frequency', 'elasticModulus', 'temperature', or 'classicCompression'.")
        else:
            csvreader = csv.reader(filePath,  delimiter=';')
            #Read header. There is just one line
            header = []
            header.append(next(csvreader))

            #Read rest of lines and buid the variable rows 
            #rows will be a nested list, with each element being one list containing the values of a line
            rows = []
            for row in csvreader:
                rowToBeAppended=[]
                for number in row:
                    try:
                        rowToBeAppended.append(float(number))
                    except ValueError:
                        rowToBeAppended.append(0)
                rows.append(rowToBeAppended)
            
            #Build object attributes, depending on the type of data being read
            if typeOfData == "frequency":
                self.testTimeSeries_elasticModulus = []
                self.testNaturalFrequencySeries = []
                for row in rows:
                    self.testTimeSeries_elasticModulus.append(row[0])
                    self.testNaturalFrequencySeries.append(row[1])
            elif typeOfData == "elasticModulus":
                self.testTimeSeries_elasticModulus = []
                self.testElasticModulusSeries = []
                for row in rows:
                    self.testTimeSeries_elasticModulus.append(row[0])
                    self.testElasticModulusSeries.append(row[1])
            elif typeOfData == "temperature":
                self.testTimeSeries_temperature = []
                self.testTemperatureSeries = []
                for row in rows:
                    self.testTimeSeries_temperature.append(row[0])
                    self.testTemperatureSeries.append(row[1])
            elif typeOfData == "classicCompression":
                self.classicCompression_Time = []
                self.classicCompression_Modulus = []
                for row in rows:
                    self.classicCompression_Time.append(row[0])
                    self.classicCompression_Modulus.append(row[1])
        
    def downsampleDataSeries(self, typeOfData, pickEveryX, firstN=0, pickEveryX_firstN=1):
        """Reduce size of a stored data series

        This method allows to downsample a data series that is too large.
        Sometimes the experimental data is too dense (too many measurement sessions) and we may downsample it to increase efficiency
        It will downsample either frequency, elastic modulus or temperature series already store in the object.

        Parameters
        ----------
        - typeOfData: str
            The type of data to be read
            Options: "frequency", "elasticModulus", "temperature"
        - pickEveryX: int
            Defines the downsampling frequency (e.g., if =10, every 10th sample will be held)
        - firstN: int
            If we want the first N samples to be sampled in a different downsampling frequency
        - pickEveryX_firstN: int
            The downsampling frequency of the first N, if defined.
        """

        if typeOfData not in [ "frequency", "elasticModulus", "temperature"]:
            print("Type of data not compatible with this method. Specify either a 'frequency', 'elasticModulus', or 'temperature'.")
        else:
            #Below we select the appropiate series to work with
            if typeOfData == "frequency":
                timeSeries= self.testTimeSeries_elasticModulus
                dataSeries = self.testNaturalFrequencySeries 
            elif typeOfData == "elasticModulus":
                timeSeries = self.testTimeSeries_elasticModulus
                dataSeries = self.testElasticModulusSeries
            elif typeOfData == "temperature":
                timeSeries = self.testTimeSeries_temperature
                dataSeries = self.testTemperatureSeries           

            #In the below loop, we will remove all elements accordingly to the downsampling rules informed in the method
            for iteration, value in enumerate(timeSeries):
                if iteration<firstN:
                    if iteration%pickEveryX_firstN!=0:
                        timeSeries.pop(iteration)
                        dataSeries.pop(iteration)
                elif iteration%pickEveryX!=0:
                    timeSeries.pop(iteration)
                    dataSeries.pop(iteration)

    def fitEvolutionCurve(self, numberOfComponents=1, verbose=True):
        """Obtain parameters of a model fitted to experimental elastic modulus measurements.

        A method to fit a model to the experimental estimates of E-modulus.
        The model is an exponential three-parameter curve.
        Multiple components comprised of three-parameters are allowed.
        It populates the attribute self.fittedEquationParameters

        Parameters
        ----------
        numberOfComponents: int, optional
            Defines the number of components of the exponential three-parameter curve
            By default, a single component is used (classical use of this model)
        verbose: boolean, optional
            Defines if textual information will be shown while executing the method.
            By default, verbose is on.
        """
        from scipy.optimize import curve_fit

        #Instantiate variable that will store fitting parameters
        self.fittedEquationParameters = {}
        for terms in np.arange(0, numberOfComponents):
            self.fittedEquationParameters.update({str(terms)+"_alfa": np.nan, str(terms)+"_beta": np.nan, str(terms)+"_tau": np.nan})

        initialGuesses = [1 for term in np.arange(0, 3*numberOfComponents)]  # Replace with value (an array)
        parametersList = [0 for term in np.arange(0, 3*numberOfComponents)]

        #Define curve to be fitted
        def func(t, *parametersList):
            #Unroll parameters list
            index = 0
            f=0
            for term in np.arange(0,len(parametersList)/3):
                α=parametersList[index]
                β=parametersList[index+1]
                𝜏=parametersList[index+2]
                f = f + α*np.exp(-(𝜏/t)**β)
                index = index + 3 
            return f
        #func = lambda t, α1, 𝜏1, β1, α2, 𝜏2, β2: (α1)*np.exp(-(𝜏1/t)**β1)+(α2)*np.exp(-(𝜏2/t)**β2)
        #Fitting procedure
        #Make the fitting
        popt, pcov = curve_fit(func, self.testTimeSeries_elasticModulus, self.testElasticModulusSeries, p0=[initialGuesses], maxfev=6000000)
        #Store fitting parameters, unrolling popt
        index = 0
        for terms in np.arange(0, numberOfComponents):
            self.fittedEquationParameters[str(terms)+"_alfa"]=popt[index]
            self.fittedEquationParameters[str(terms)+"_beta"]=popt[index+1]
            self.fittedEquationParameters[str(terms)+"_tau"]=popt[index+2]
            index =+ 3
        
        self.fittedStandardError=np.sqrt(np.diag(pcov))

    def computeFromFittedModel(self, timeSeries=None, whatToCompute="all", methodForDerivative="analytical", callContext="general"):
        """Compute information from the fitted elastic modulus model.
        
        Use the fitted model to derive elastic modulus, advancement degree and derivative of advancement degree curve

        Parameters
        ----------
        timeSeries: list of int, optional
            The time series which will be used to perform those computations
            If not given, the experimental time series associated to elastic modulus will be used
        whatToCompute: str, optional
            Parameter that specifies what information is to be computed from the model.
            Options are: "all", "elasticModulus", "advancementDegree", "derivativeAdvancementDegree".
            Default is "all", which computes all possible information
        methodForDerivative: str, optional
            Defines how the derivative will be computed.
            Options are: "analytical" or "numerical".
            Default is "analytical".
        callContext: str, optional
            Defines if the context in which this method is called.
            Options are: "general", "activationEnergyCalculation".
            Different context may allow for optimized calculations or setting of specific attributes.
            Default is "general".

        Raises
        ------
        AttributeError
            If self.fittedEquationParameters is not set previous to this method.
        """

        #Initial verifications
        if self.fittedEquationParameters is None:
            raise AttributeError("Model was not fitted. Fit it first then use this method.")
        
        #Computations of ultimate value from the fitted model's parameter 
        ultimateValue=0
        for value in self.fittedEquationParameters:
            if "alfa" in value:
                ultimateValue = ultimateValue + self.fittedEquationParameters[value]

        #Initialize some variables
        #If it is a general call, perform the usual computations
        if callContext == "general":
            #Initial validation of input parameters
            if timeSeries is None:
                self.fittedTimeSeries = np.array(self.testTimeSeries_elasticModulus)
            else:
                self.fittedTimeSeries = np.array(timeSeries)
            #Check what needs to be computed and initialize variables
            if whatToCompute == "all" or whatToCompute == "elasticModulus":
                self.fittedElasticModulus=np.zeros(len(self.fittedTimeSeries))
            if whatToCompute == "all" or whatToCompute == "advancementDegree":
                self.fittedAdvancementDegree=np.zeros(len(self.fittedTimeSeries))
            if whatToCompute == "all" or whatToCompute == "derivativeAdvancementDegree":
                self.fittedDerivativeAdvancementDegree=np.zeros(len(self.fittedTimeSeries))
            #Do the actual computation of required information
            for terms in np.arange(0, len(self.fittedEquationParameters)/3):
                α=self.fittedEquationParameters[str(int(terms))+"_alfa"]
                β=self.fittedEquationParameters[str(int(terms))+"_beta"]
                𝜏=self.fittedEquationParameters[str(int(terms))+"_tau"]
                if whatToCompute == "all" or whatToCompute == "elasticModulus":
                    self.fittedElasticModulus = self.fittedElasticModulus + α*np.exp(-(𝜏/self.fittedTimeSeries)**β)
                if whatToCompute == "all" or whatToCompute == "advancementDegree":
                    self.fittedAdvancementDegree = self.fittedAdvancementDegree + (α/ultimateValue)*np.exp(-(𝜏/self.fittedTimeSeries)**β)
                if whatToCompute == "all" or whatToCompute == "derivativeAdvancementDegree":
                    if methodForDerivative == "analytical":
                        self.fittedDerivativeAdvancementDegree = self.fittedDerivativeAdvancementDegree + (α/ultimateValue)*(β*(𝜏**β))*(1/(self.fittedTimeSeries**(β+1)))*np.exp(-(𝜏/self.fittedTimeSeries)**β)
            if methodForDerivative == "numerical":
                self.fittedDerivativeAdvancementDegree = np.gradient(self.fittedAdvancementDegree,self.fittedTimeSeries)
            elif methodForDerivative != "analytical":
                raise Exception("The methodForDerivative informed is not valid. It should be either 'analytical' or 'numerical'")
        elif callContext == "activationEnergyCalculation":
            #If it is a call within the Speed or Derivative Methods, perform its custom computations
            #In this context, we can only be interested in computing the derivative of advancement degree

            #As opposed to the general call, in this context we don't set self.fittedTimeSeries
            #We substitute that for a "shadow copy" for the speed method: self.fittedTimeSeries_EaCalculation
            if timeSeries is None:
                self.fittedTimeSeries_EaCalculation = np.array(self.testTimeSeries_elasticModulus)
            else:
                self.fittedTimeSeries_EaCalculation = np.array(timeSeries)
            self.derivativeAdvancementDegree_EaCalculation = np.zeros(len(self.fittedTimeSeries_EaCalculation))
            if methodForDerivative == "analytical":
                for terms in np.arange(0, len(self.fittedEquationParameters)/3):
                    α=self.fittedEquationParameters[str(int(terms))+"_alfa"]
                    β=self.fittedEquationParameters[str(int(terms))+"_beta"]
                    𝜏=self.fittedEquationParameters[str(int(terms))+"_tau"]
                    self.derivativeAdvancementDegree_EaCalculation = self.derivativeAdvancementDegree_EaCalculation + (α/ultimateValue)*(β*(𝜏**β))*(1/(self.fittedTimeSeries_EaCalculation**(β+1)))*np.exp(-(𝜏/self.fittedTimeSeries_EaCalculation)**β)
            elif methodForDerivative == "numerical":
                self.derivativeAdvancementDegree_EaCalculation = np.gradient(self.fittedAdvancementDegree,self.fittedTimeSeries_EaCalculation)
            else:
                raise Exception("The methodForDerivative informed is not valid. It should be either 'analytical' or 'numerical'")
        else:
            raise Exception("Parameter isSecondaryInActivationEnergy_SpeedMethod is invalid. Should be boolean.")

        #Identifies which time series was used to compute the fitted attributes and store this information
        if np.array(timeSeries==self.testTimeSeries_elasticModulus).all():
            #Then, experimental time series associated to elastic modulus was used
            self.typeOfTimeSeriesUsedInFitting="experimental_elasticModulus"
        elif np.array(timeSeries==self.testTimeSeries_temperature).all():
            #Then, experimental time series associated to temperature measurement was used
            self.typeOfTimeSeriesUsedInFitting="experimental_temperature"
        else:
            #Another time series not configure (probably arbitraly, was used)
            self.typeOfTimeSeriesUsedInFitting="custom"

    def computeEquivalentAge(self, activationEnergySeries, Tref, thresholdDegree=0.01):
        """Compute the equivalent age of the self.fittedTimeSeries in a reference temperature Tref.
        
        Using a given activation energy series, compute what is the equivalent age of each time instant contained in self.fittedTimeSeries attribute.
        Since the activation energy value vary according to the degree of advancement,
        we need to make sure the activation energy, temperature and time steps in the equivalent age calculation all refer to the same degree of advancement.
        
        For that, we assess the attribute self.fittedAdvancementDegree, which contains the correspondence between time instants and degrees.
        For a given specimen, the degrees containing in self.fittedAdvancementDegree may not coincide with the degrees in the activation energy series.
        For that, we use a linear interpolation to get an approximate activation energy and temperature for each degree and time instant.

        Parameters
        ----------
        activationEnergySeries: 2D Numpy array
            Series of degree of advancement VS activation energy values.
            The first column should contain the degree of advancement values.
            The second column should contain the activation energy values.
        Tref: float
            A reference temperature in relation to which the equivalent age is computed.
        thresholdDegree: float, optional
            The minimum degree of advancement below which all values from the activation energy series will be ignored.
            This is necessary because for very low degrees (usually <1%), the test provides unrealistic values
            that, in equivalent age calculation, are related to large propagations of error. 
        
        Raises
        ------
        """
        import copy
        #Make a deep copy so we don't edit activationEnergySeries which may be an external variable in our code
        activationEnergySeriesCopy = copy.deepcopy(activationEnergySeries)

        '''
        #Original implementation
        #Pre-process the activation energy series to remove any spuriour values that may happen
        #in very early hydration degrees, in which thermal imbalance and other effects may lead
        #to clearly wrong activation energy values
        indicesToProcess = [index for index, value in enumerate(activationEnergySeriesCopy) if value[0]<thresholdDegree]
        for index in indicesToProcess:
            #activationEnergySeriesCopy[index,1]=activationEnergySeriesCopy[max(indicesToProcess)+1,1]
            #Now we disregard all entries that are associated to degrees smaller than threshold
            activationEnergySeriesCopy=activationEnergySeriesCopy[max(indicesToProcess)+1:,:]

        #Initialize the attribute self.equivalentAge
        self.equivalentAge = np.zeros(len(self.fittedTimeSeries))
        t_previous = 0
        for i, t in enumerate(self.fittedTimeSeries):
            #First we build the maps that will give, for each instant, the associated Eact and temperature
            Eac=np.interp(self.fittedAdvancementDegree[i], np.concatenate((np.array([0]),activationEnergySeriesCopy[:,0])), np.concatenate((np.array([activationEnergySeriesCopy[0,1]]),activationEnergySeriesCopy[:,1])))
            Ti=np.interp(t, np.concatenate((np.array([0]),self.testTimeSeries_temperature)), np.concatenate((np.array([self.testTemperatureSeries[0]]),self.testTemperatureSeries)))
            #Compute the time step Δt associated to the current instant t
            Δt = t - t_previous
            t_previous = t
            #Compute the equivalent age associated to the current self.fittedTimeSeries instant
            if i == 0:
                self.equivalentAge[i]=np.exp((-Eac/R)*((1/(Ti+273.15))-(1/(Tref+273.15))))*Δt
            else:
                self.equivalentAge[i]=self.equivalentAge[i-1]+np.exp((-Eac/R)*((1/(Ti+273.15))-(1/(Tref+273.15))))*Δt     
        '''
        #Second attempt
        #Pre-process the activation energy series to remove any spuriour values that may happen
        #in very early hydration degrees, in which thermal imbalance and other effects may lead
        #to clearly wrong activation energy values
        indicesToDisregard = [index for index, value in enumerate(activationEnergySeriesCopy) if value[0]<thresholdDegree]
        #Now we disregard all entries that are associated to degrees smaller than threshold
        #The below try-except handles the case in which no index is to be disregarded (len(indicesToDisregard)=0)
        try:
            activationEnergySeriesCopy=activationEnergySeriesCopy[max(indicesToDisregard)+1:,:]
        except ValueError:
            pass

        #Now, check inside self.fittedTimeSeries and self.fittedAdvancementDegree, which entries are associated
        #to degrees smaller than the thresholdDegree
        indicesToDisregard = [index for index, value in enumerate(self.fittedAdvancementDegree) if value<thresholdDegree]
        #The below try-except handles the case in which no index is to be disregarded (len(indicesToDisregard)=0)
        try:
            startIndex = max(indicesToDisregard)+1
        except ValueError:
            startIndex = 0
        #Initialize the attribute self.equivalentAge, with the correct size
        self.equivalentAge = np.zeros(len(self.fittedTimeSeries[startIndex:]))
        t_previous = 0
        #Only loop through values in self.fittedTimeSeries not included in indicesToDisregard!
        for i, t in enumerate(self.fittedTimeSeries[startIndex:]):
            #First we build the maps that will give, for each instant, the associated Eact and temperature
            #Here we need to correct the indices in self.fittedAdvancementDegree
            Eac=np.interp(self.fittedAdvancementDegree[i+startIndex], np.concatenate((np.array([0]),activationEnergySeriesCopy[:,0])), np.concatenate((np.array([activationEnergySeriesCopy[0,1]]),activationEnergySeriesCopy[:,1])))
            Ti=np.interp(t, np.concatenate((np.array([0]),self.testTimeSeries_temperature)), np.concatenate((np.array([self.testTemperatureSeries[0]]),self.testTemperatureSeries)))
            #Compute the time step Δt associated to the current instant t
            Δt = t - t_previous
            t_previous = t
            #Compute the equivalent age associated to the current self.fittedTimeSeries instant
            if i == 0:
                self.equivalentAge[i]=np.exp((-Eac/R)*((1/(Ti+273.15))-(1/(Tref+273.15))))*Δt
            else:
                self.equivalentAge[i]=self.equivalentAge[i-1]+np.exp((-Eac/R)*((1/(Ti+273.15))-(1/(Tref+273.15))))*Δt    

    def interpolateTemperatureSeries(self, interpolatingTimeSeries):
        '''This function interpolates the experimental temperature series to a specific time series.
        
        Parameters
        ----------
        interpolatingTimeSeries: 1D Numpy array
            Time series for which we want interpolated values of the temperature series
        
        Returns
        ------
        interpolatedTemperatureSeries: 1D Numpy array
            Temperature series interpolated at the desired time series.
        '''
        
        #Initialize the attributes with the correct length
        interpolatedTemperatureSeries=np.zeros(len(interpolatingTimeSeries)) 
        for index, t in enumerate(interpolatingTimeSeries):
            #We need to append to the beggining of the series a time=0 and a corresponding temperature
            #Otherwise if interpolating time t lesser than the smallest time in self.testTimeSeries_temperature
            #is being interpolated, it might return weird results (not clear if np.interp can handle that automatically)
            #So we just append to self.testTimeSeries_temperature a time = 0 
            #and to self.testTemperatureSeries its first temperature retrieved from self.testTemperatureSeries[0]
            interpolatedTemperatureSeries[index]=np.interp(t, np.concatenate((np.array([0]),self.testTimeSeries_temperature)), np.concatenate((np.array([self.testTemperatureSeries[0]]),self.testTemperatureSeries)))

        return interpolatedTemperatureSeries

class experiment_EMMARM ():
    """
    A class that defines all a EMMARM test specimen

    Attributes:
    - testId: str
        String that identifies the test
    - testSpecimensList: list, of EMMARMspecimen objects
        Contains the EMMARM specimens that are part of a single experiment
    - activationEnergy_DerivativeMethod: dict, of 3D lists
        Activation energy series will be contained in a dictionary whose index will be formatted as "primarySpecimenID". 
        Each activation energy entry will be associated to a 5D list.
        The general formating of a single line of a 5D list associated to a given dictionary of this list will be:
            - {"primarySpecimenID": np.array([degree, Eact, C, [[x1,x2,x3,...,xn],[y1,y2,y3,...,yn]],r2])}
            - The first column contains the advancement degree values, the second column contains the respective activation energy estimates, and the third column contains the linear constant values.
            - The fourth colum will store a 2D list of point values (x,y) associated to each "n" specimens that compose the Arrhenius plot of each degree of advancement, so the Arrhenius plots can be reconstructed if needed. These lists will have dimension of "n" (as each pair of (x,y) is associated to a single specimen)
            - The fifth column will store the r2 (correlation coefficient) parameter associated to the fitting of the current line
    - activationEnergy_SpeedMethod: dict, of 2D Numpy arrays
        Activation energy series will be contained in a dictionary whose index will be formatted as "PrimarySpecimenID"-"SecondarySpecimenID". 
        Each activation energy entry will be associated to a 2D Numpy array, in which the first column contains the advancement degree values and the second column contains the respective activation energy estimates.
    - self.activationEnergy_SpeedMethod_averaged: 2D Numpy array
        Contains an averaged activation energy series based on all estimate pairs obtained from the Speed Method.
        To do the averaging, multiple series need to be interpolated to same activation energy degree series.
        The first column of this 2D Numpy array contains the advancement degree and the second contains the activation energy values.

    Methods:
    - ...
    """

    def __init__(self, testId, testSpecimensList):
        """Initialization of the class.

        A experiment is define by one or more specimens.
        They may be replicates tested at the same time or sequentially, either at single or multiple temperatures.

        Parameters
        ----------
        - testId: str
            Identification of this test
        - specimenList: list, of EMMARMspecimen objects
            Contains the EMMARM specimens that are part of a single experiment
        """
        self.testId=testId
        self.testSpecimensList=testSpecimensList

        #Attributes set by methods
        self.activationEnergy_DerivativeMethod = None 
        self.activationEnergy_SpeedMethod = None
        self.activationEnergy_SpeedMethod_averaged = None
        
    #Methods related to the computation of activation energy via Derivative of Speed Method
    def initializeActivationEnergy_DerivativeMethod(self):
        """Initialize the attribute self.activationEnergy_DerivativeMethod

        All possible combinations based on the test specimens comprising the test are calculated, and all the respective activation energy series are initialized with zero values and the correct length.
        For "n" specimens, there will be n combinations (because we are looking for the total number of combinations with a unique Primary series - the order of the Secondary series does not matter).
        This will be the length of the self.activationEnergy_DerivativeMethod attribute.
        """
        numberOfSpecimens=len(self.testSpecimensList)
        self.activationEnergy_DerivativeMethod={}
        for indexPrimary in np.arange(0, numberOfSpecimens):
            #Select the specimen that will be the Primary series in the current combination
            idPrimary=self.testSpecimensList[indexPrimary].specimenID
            #Define an object called generalLine which will store the general format of a single line
            generalLine = [0,0,0,[[0 for specimen in self.testSpecimensList],[0 for specimen in self.testSpecimensList]],0]
            #Initialize the activation energy variable for this combination
            self.activationEnergy_DerivativeMethod.update({idPrimary: [generalLine for degree in self.testSpecimensList[indexPrimary].fittedAdvancementDegree]})

    def computeActivationEnergy_DerivativeMethod(self):
        ##Initialize the attribute to store the activation energy estimates
        #Define all possible combinations of the tests
        #In the current implementation, we may only allow to use all tests at the same time
        #But we could also implement a way to use only part of them (not so interesting imo)
        #For "n" tests, there would be "n" combinations: each time, one test is taken as Primary and the others as Secondary
        #{"primarySpecimenID": np.array([degree, Eact, C, [[x1,x2,x3,...,xn],[y1,y2,y3,...,yn]])}
        self.initializeActivationEnergy_DerivativeMethod()

        numberOfSpecimens=len(self.testSpecimensList)
        for indexPrimary in np.arange(0, numberOfSpecimens):
            primarySpecimen=self.testSpecimensList[indexPrimary]
            #Set the appropriate temperature and derivative series
            #Set Primary series temperature series
            if primarySpecimen.typeOfTimeSeriesUsedInFitting == "experimental_temperature":
                primarySpecimen.temperatureSeries_ForEaCalculation=primarySpecimen.testTemperatureSeries
            else:
                primarySpecimen.temperatureSeries_ForEaCalculation = primarySpecimen.interpolateTemperatureSeries(primarySpecimen.fittedTimeSeries)
            primarySpecimen.derivativeAdvancementDegree_EaCalculation=primarySpecimen.fittedDerivativeAdvancementDegree
            #Interpolate Secondary series temperature and derivative degree of advancement to Primary series degrees of advancement
            for indexSecondary in np.arange(0,numberOfSpecimens):
                #We just want to select Secondary series, not the Primary, so indexSecondary != indexPrimary
                #if indexSecondary != indexPrimary:
                secondarySpecimen=self.testSpecimensList[indexSecondary]
                self.interpolateSecondarySeries_EaCalculation(primarySpecimen, secondarySpecimen)
            #Define the fitting curve
            from scipy.optimize import curve_fit
            from scipy.stats import linregress
            def arrheniusCurve(T, Eact, C):
                return C - (Eact/R)*(1/(T+273.15))
            #For all degree of advancement, fit the arrheniusCurve to the set of (x,y) according to the Derivative Method
            for index, degree in enumerate(primarySpecimen.fittedAdvancementDegree):
                '''
                #First attempted alternative: using curve_fit
                #Build lists with fitting points
                x = [specimen.temperatureSeries_ForEaCalculation[index] for specimen in self.testSpecimensList]
                y = [np.log(specimen.derivativeAdvancementDegree_EaCalculation[index]) for specimen in self.testSpecimensList]
                self.activationEnergy_DerivativeMethod[primarySpecimen.specimenID][index]=[degree, fittedEact, fittedC, [x,y]]
                '''
                #Second alternative: using linregress. Preferred due to simplicity and speed
                x = [1/(specimen.temperatureSeries_ForEaCalculation[index]+273.15) for specimen in self.testSpecimensList]
                y = [np.log(specimen.derivativeAdvancementDegree_EaCalculation[index]) for specimen in self.testSpecimensList]
                slope, intercept, r, p, se = linregress(x, y)
                self.activationEnergy_DerivativeMethod[primarySpecimen.specimenID][index]=[degree, -slope*R, intercept, [x,y],r]

    def getActivationEnergySeries_DerivativeMethod(self, dictKeyword):
        """A method to retrieve a 2D np.array containing the activation energy series.
        
        Parameters
        ----------
        dictKeyword: str
            A string that identifies the desired activation energy series.
            It is the ID of the specimen chosen as Primary series for the desired calculation.

        Returns
        -------
        activationEnergySeries: 2D ndarray
            A 2D Numpy array with first column containing degrees of advancement and second column containing activation energy values.
        """ 
        activationEnergySeries=np.zeros((len(self.activationEnergy_DerivativeMethod[dictKeyword]),2))
        for index, line in enumerate(self.activationEnergy_DerivativeMethod[dictKeyword]):
            activationEnergySeries[index,0]=line[0]
            activationEnergySeries[index,1]=line[1]
        
        return activationEnergySeries

    #Methods related to the computation of activation energy via Speed Method
    def initializeActivationEnergy_SpeedMethod(self):
        """Initialize the attribute self.activationEnergy_SpeedMethod

        All possible combinations based on the test specimens comprising the test are calculated, and all the respective activation energy series are initialized with zero values and the correct length.
        For "n" specimens, there will be n!/(n-2)! combinations (because we are looking for the total number of ordered pairs).
        This will be the length of the self.activationEnergy_SpeedMethod attribute.
        
        """
        numberOfSpecimens=len(self.testSpecimensList)
        self.activationEnergy_SpeedMethod={}
        for indexPrimary in np.arange(0, numberOfSpecimens):
            for indexSecondary in np.arange(indexPrimary+1, numberOfSpecimens):
                #Select the specimens that will be combined in a pair for Speed Method
                idPrimary=self.testSpecimensList[indexPrimary].specimenID
                idSecondary=self.testSpecimensList[indexSecondary].specimenID
                #Compute the length of the activation energy to be calculated
                #It will be equal to length of Primary series duration, as it will interpolate Secondary series
                lengthOfAdvancementDegreeSeries=len(self.testSpecimensList[indexPrimary].fittedAdvancementDegree)
                #Initialize the activation energy variable for this pair
                self.activationEnergy_SpeedMethod.update({idPrimary+"-"+idSecondary: np.zeros((lengthOfAdvancementDegreeSeries,2))})
                #Also, do the same for the inverse pair (i.e., Secondary-Primary)
                lengthOfAdvancementDegreeSeries=len(self.testSpecimensList[indexSecondary].fittedAdvancementDegree)
                self.activationEnergy_SpeedMethod.update({idSecondary+"-"+idPrimary: np.zeros((lengthOfAdvancementDegreeSeries,2))})

    def computeActivationEnergy_SpeedMethod(self, temperatureSeriesToUse="auto"):
        """Compute apparent activation energy via Speed Method

        Compute the activation energy for pairs of isothermal EMM-ARM tests that comprise the current experiment.
        All possible ordered pairs are evaluated, so both combinations specimenA-specimenB and specimenB-specimenA are made.
        The order matters because the first series will always be taken as the reference.
        The first series will be referenced as Primary.
        The second series (referenced as Secondary) will temperature and advacement degree series will be interpolated to fit the hydration degrees of the Primary series.

        #TODO: Perhaps this contraint may be enough to derive a sub-class for Speed Method from experiment_EMMARM class, using deep copies
        
        Parameters:
            temperatureSeriesToUse: str, optional
                Defines which temperature series to use in the method for the Primary specimen.
                Options are "interpolated" or "experimental"
                Using "interpolated" will use the interpolated temperature series in self.interpTemperatureSeries
                Using "experimental" will use the experimental temperature series in self.testTemperatureSeries

        """
        #Initialize the attribute self.initializeActivationEnergy_SpeedMethod
        self.initializeActivationEnergy_SpeedMethod()

        #Compute activation energy for each pair of specimens
        numberOfSpecimens=len(self.testSpecimensList)
        for indexPrimary in np.arange(0, numberOfSpecimens):
            for indexSecondary in np.arange(indexPrimary+1, numberOfSpecimens):
                #First, we select the specimens
                primarySpecimen=self.testSpecimensList[indexPrimary]
                secondarySpecimen=self.testSpecimensList[indexSecondary]
                #Select the appropriate temperature series of the Primary series
                #This is based on what was used to compute its own parameters 
                #from the fitted curve when using the method computeFromFittedModel()
                if primarySpecimen.typeOfTimeSeriesUsedInFitting == "experimental_temperature":
                    primarySpecimen.temperatureSeries_ForEaCalculation=primarySpecimen.testTemperatureSeries
                else:
                    primarySpecimen.temperatureSeries_ForEaCalculation = primarySpecimen.interpolateTemperatureSeries(primarySpecimen.fittedTimeSeries)
                #Then, we interpolated the Secondary series derivative of advancement degree and temperature
                self.interpolateSecondarySeries_EaCalculation(primarySpecimen,secondarySpecimen)
                #Now we compute the apparent activation energy via Speed Method
                idPrimary=primarySpecimen.specimenID
                idSecondary=secondarySpecimen.specimenID
                activationEnergyIndex=idPrimary+"-"+idSecondary
                for index, advancentDegree in enumerate(primarySpecimen.fittedAdvancementDegree):
                    self.activationEnergy_SpeedMethod[activationEnergyIndex][index, 0]=advancentDegree
                    primT = primarySpecimen.temperatureSeries_ForEaCalculation[index]
                    secdT = secondarySpecimen.temperatureSeries_ForEaCalculation[index]
                    primDervAlfa = primarySpecimen.fittedDerivativeAdvancementDegree[index]
                    secdTDervAlfa = secondarySpecimen.derivativeAdvancementDegree_EaCalculation[index]
                    self.activationEnergy_SpeedMethod[activationEnergyIndex][index, 1]=((-R/((1/(primT+273.15))-(1/(secdT+273.15))))*(np.log(primDervAlfa)-np.log(secdTDervAlfa)))
                #Now, we do the same but with the inverse pair of specimens
                primarySpecimen=self.testSpecimensList[indexSecondary]
                secondarySpecimen=self.testSpecimensList[indexPrimary]
                self.interpolateSecondarySeries_EaCalculation(primarySpecimen,secondarySpecimen)
                idPrimary=primarySpecimen.specimenID
                idSecondary=secondarySpecimen.specimenID
                activationEnergyIndex=idPrimary+"-"+idSecondary
                for index, advancentDegree in enumerate(primarySpecimen.fittedAdvancementDegree):
                    self.activationEnergy_SpeedMethod[activationEnergyIndex][index, 0]=advancentDegree
                    primT = primarySpecimen.testTemperatureSeries[index]
                    secdT = secondarySpecimen.temperatureSeries_ForEaCalculation[index]
                    primDervAlfa = primarySpecimen.fittedDerivativeAdvancementDegree[index]
                    secdTDervAlfa = secondarySpecimen.derivativeAdvancementDegree_EaCalculation[index]
                    self.activationEnergy_SpeedMethod[activationEnergyIndex][index, 1]=((-R/((1/(primT+273.15))-(1/(secdT+273.15))))*(np.log(primDervAlfa)-np.log(secdTDervAlfa)))
        #Compute the average of all activation energy pairs
        self.averageActivationEnergy_SpeedMethod()    
        
        '''
        #TODO: Check if this is really necessary. If it is, implement after activation energy calculation
        #Ignore spurious values originated from numerical errors
        activationEnergyValidBounds=[-10, 200]
        advancementDegreeValidBounds=[1e-60,0.99]
        indicesToExclude=[]
        for iteration, degree in enumerate(hydrationDegree[0]):
            if degree<hydrationDegreeValidBounds[0] or degree>hydrationDegreeValidBounds[1]:
                indicesToExclude.append(iteration)
        activationEnergy=np.delete(activationEnergy, indicesToExclude)
        hydrationDegree[0]=np.delete(hydrationDegree[0], indicesToExclude)
        indicesToExclude=[]
        for iteration, value in enumerate(activationEnergy):
            if value<activationEnergyValidBounds[0] or value>activationEnergyValidBounds[1]:
                indicesToExclude.append(iteration)
        activationEnergy=np.delete(activationEnergy, indicesToExclude)
        hydrationDegree[0]=np.delete(hydrationDegree[0], indicesToExclude)
        '''

    def averageActivationEnergy_SpeedMethod(self):
        '''Method to average the collection of activation energy estimates of Speed Method
        
        Speed method provides various activation energy estimates because it only works with a pair of isothermal test. This method will average them all. To do that, we must first interpolate some of them to the same degree of advancement values, such that averaging can be done
        '''

        from scipy import interpolate
        #To facilitate computation, store the advancement degree VS activation energy series, originally in a dictionary, in separate collection of lists of α VS Eact
        α = [[] for series in self.activationEnergy_SpeedMethod]
        Eact = [[] for series in self.activationEnergy_SpeedMethod]
        for iteration, series in enumerate(self.activationEnergy_SpeedMethod):
            for line in self.activationEnergy_SpeedMethod[series]:
                α[iteration].append(line[0])
                Eact[iteration].append(line[1])
        #Discover the limits of interpolation based on minimal and maximum values of the original Eact
        max_min_α=max([min(series) for series in α])
        min_max_α=min([max(series) for series in α])
        #Define the new degree of advancement series
        new_α = np.arange(max_min_α,min_max_α,0.005)
        interpolated_Eact = [np.zeros(len(new_α)) for series in self.activationEnergy_SpeedMethod]
        self.activationEnergy_SpeedMethod_averaged = np.zeros((len(new_α),2))
        #Now interpolate original Eact to same advancement degrees so we can average them
        for iteration, series in enumerate(self.activationEnergy_SpeedMethod):
            interpolationFunction = interpolate.interp1d(α[iteration], Eact[iteration])
            interpolated_Eact[iteration] = interpolationFunction(new_α)
        #Average all three series
        for iteration, degree in enumerate(new_α):
            self.activationEnergy_SpeedMethod_averaged[iteration,0]=degree
        for series in interpolated_Eact:
            for iteration, element in enumerate(series):
                self.activationEnergy_SpeedMethod_averaged[iteration,1]=self.activationEnergy_SpeedMethod_averaged[iteration, 1]+element
        self.activationEnergy_SpeedMethod_averaged[:,1]=self.activationEnergy_SpeedMethod_averaged[:,1]/len(interpolated_Eact)

    def getActivationEnergySeries_SpeedMethod(self, dictKeyword):
        """A method to retrieve a 2D np.array containing the activation energy series.
        
        Parameters
        ----------
        dictKeyword: str
            A string that identifies the desired activation energy series.
            It is the ID of the specimen chosen as Primary series for the desired calculation
            followed by the ID of the specimen chosen as Secondary series.

        Returns
        -------
        activationEnergySeries: 2D ndarray
            A 2D Numpy array with first column containing degrees of advancement and second column containing activation energy values.
        """ 
        activationEnergySeries=np.zeros((len(self.activationEnergy_SpeedMethod[dictKeyword]),2))
        for index, line in enumerate(self.activationEnergy_SpeedMethod[dictKeyword]):
            activationEnergySeries[index,0]=line[0]
            activationEnergySeries[index,1]=line[1]
        
        return activationEnergySeries

    #General methods used in both activation energy methods
    def interpolateSecondarySeries_EaCalculation(self, primarySpecimen, secondarySpecimen, methodForInterpolation = 'approximate'):
        '''This method interpolates the derivative of advancement degree and temperature of the Secondary series 
        
        In the context of the Speed Method for calculating apparent activation energy, it is necessary to interpolate the derivatives of advancement degree and temperature series of the Secondary series so they refer to the same hydration degree instants of the Primary series.
        This method does this and store the results in the proper attributes of the Secondary series object.
        These attributes are temporary for a given specimen, such that they will be written over every time this method is called with such specimen. 
        So they are only updated to the last call of this method.

        Parameters:
            primarySpecimen: EMMARMspecimen object
                The Primary specimen in the context of the Speed Method
            secondarySpecimen: EMMARMspecimen object
                The Secondary specimen in the context of the Speed Method
            methodForInterpolation: str, optional
                Defines the method for interpolation.
                Options are:
                    "approximate" - values are interpolated from linear interpolation
                    "implict" - values are interpolated from implicit solving the fitted model equation
                Default is "approximate" which is quicker

        '''
        
        interpolatedTimes=[]
        #Compute the times of Secondary series in which its hydration degree values are the sames as the Primary series
        if methodForInterpolation=="approximate":
            #This method is a simple linear interpolation so we know in which times of the Secondary test we achieved the same hydration degrees as those achieved in the measurement times of the Primary test
            interpolatedTimes = np.interp(primarySpecimen.fittedAdvancementDegree,secondarySpecimen.fittedAdvancementDegree, secondarySpecimen.fittedTimeSeries)
        elif methodForInterpolation=="implicit":
            from scipy import optimize
            for iteration, α in enumerate(primarySpecimen.fittedAdvancementDegre): 
                #Define the function f which is our objective function
                def f(x):
                    '''
                    This function computes the difference between the current advancement degree α and the degree in the Secondary series associated to an instant x. 
                    The value x that minimizes f will be the one we will store in our interpolatedTimes vector
                    '''
                    functionValue=0
                    for value in secondarySpecimen.fittedEquationParameters:
                        if "alfa" in value:
                            ultimateValue = ultimateValue + secondarySpecimen.fittedEquationParameters[value]
                    for terms in np.arange(0, len(secondarySpecimen.fittedEquationParameters)/3):
                        α=secondarySpecimen.fittedEquationParameters[str(int(terms))+"_alfa"]
                        β=secondarySpecimen.fittedEquationParameters[str(int(terms))+"_beta"]
                        𝜏=secondarySpecimen.fittedEquationParameters[str(int(terms))+"_tau"]
                        functionValue = functionValue + (α/ultimateValue)*np.exp(-(𝜏/x)**β)
                    return functionValue
                #Define lower and upper bounds for optimization search
                if iteration == 0:
                    lowerBound=0
                    higherBound=primarySpecimen.fittedTimeSeries[0]
                    #initialGuess = 1e-60
                else:
                    lowerBound = interpolatedTimes[iteration-1]
                    higherBound = lowerBound+1
                    #initialGuess = t[iteration-1]
                #Perform optimization (implicit solving)
                #t.append(fsolve(f, initialGuess)[0])
                #t.append(minimize_scalar(f, bounds=(lowerBound, higherBound)).x)
                interpolatedTimes.append(optimize.root_scalar(f, bracket=[lowerBound, higherBound], method='brentq').root)
        else:
            raise Exception("Error. Method for activation energy calculation invalid. Script terminated.")
        
        #Now, use the interpolated times to find interpolated derivatives of advancement degrees
        secondarySpecimen.computeFromFittedModel(timeSeries=interpolatedTimes, whatToCompute='derivativeAdvancementDegree', callContext="activationEnergyCalculation")
        #And the interpolated temperatures
        secondarySpecimen.temperatureSeries_ForEaCalculation=np.interp(interpolatedTimes, secondarySpecimen.testTimeSeries_temperature,secondarySpecimen.testTemperatureSeries)

if __name__=="__main__":
    cementPaste_Temperature_cp16 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp16_temperature.csv')
    cementPaste_Modulus_cp16 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp16_elasticmodulus.csv')
    cementPaste_Temperature_cp31 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp31_temperature.csv')
    cementPaste_Modulus_cp31 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp31_elasticmodulus.csv')
    cementPaste_Temperature_cp40 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp40_temperature.csv')
    cementPaste_Modulus_cp40 = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp40_elasticmodulus.csv')
    fileConcreteEModulus = open('C:/Users/renan/OneDrive/Renan/Engenharia Civil/Artigos/2023-ToBeDefined-ActivationEnergy/02-DATA_PROCESSING/rawData/EMMcp31_temperature.csv')

    EMMcp16 = EMMARMspecimen("EMMcp16")
    EMMcp16.readProcessedTestData(cementPaste_Modulus_cp16, "elasticModulus")
    EMMcp16.readProcessedTestData(cementPaste_Temperature_cp16, "temperature")
    EMMcp31 = EMMARMspecimen("EMMcp31")
    EMMcp31.readProcessedTestData(cementPaste_Modulus_cp31, "elasticModulus")
    EMMcp31.readProcessedTestData(cementPaste_Temperature_cp31, "temperature")
    EMMcp40 = EMMARMspecimen("EMMcp40")
    EMMcp40.readProcessedTestData(cementPaste_Modulus_cp40, "elasticModulus")
    EMMcp40.readProcessedTestData(cementPaste_Temperature_cp40, "temperature")

    EMMcp16.fitEvolutionCurve(numberOfComponents=2)   
    EMMcp16.computeFromFittedModel(timeSeries = EMMcp16.testTimeSeries_temperature) 
    EMMcp31.fitEvolutionCurve(numberOfComponents=2)
    EMMcp31.computeFromFittedModel(timeSeries = EMMcp31.testTimeSeries_temperature) 
    EMMcp40.fitEvolutionCurve(numberOfComponents=2)
    EMMcp40.computeFromFittedModel(timeSeries = EMMcp40.testTimeSeries_temperature) 

    exp1=experiment_EMMARM("test1",[EMMcp16, EMMcp31, EMMcp40])
    exp1.computeActivationEnergy_SpeedMethod()
    exp1.computeActivationEnergy_DerivativeMethod()
    #computeActivationEnergy_SpeedMethod()
    print("The End")
    