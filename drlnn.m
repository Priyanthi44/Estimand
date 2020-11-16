side=1;
num=1;
items=1000;

Data = Dataset(side,num,items);

%% create training and testing matrices

[entries, attributes] = size(Data);

entries_breakpoint = round(entries*.90); 
%set breakpoint for training and testing data at 90% of dataset

inputlayersize=2;

outputlayersize=attributes-inputlayersize;

trainingdata = Data(1:entries_breakpoint,:);
%truncate first 90% entries for training data

trainingdata_inputs = trainingdata(:,1:inputlayersize);
%90%x9 matrix input training data

trainingdata_outputs = trainingdata(:,inputlayersize+1:end); 
%90:1 matrix output training data

testingdata = Data(entries_breakpoint:end,:); 
%truncate last 10 entries for testing data

        testingdata_inputs= testingdata(:,1:inputlayersize); %10:9 matrix input testing data

        testingdata_outputs= testingdata(:,inputlayersize+1:end); %10:1 matrix output testing data