#Integrated Environmental Health and Safety index (IEHS)

This is a design tool for chemical factory designers in the process flowsheet and conceptual design phase in the larger phase of Front End Engineering and Design (FEED). It can be used throughout the stage from an Input-Output structure to the final flowsheet alongside a rigorous technoeconmic study. It runs on MatLab and Octave with minimal system requirements. This is a redesigned index from the original EHS index by Koller et al. (2000). The index is broken up into 13 different scenarios, called parameters, that are tested for hazards. These are queried at different scales in a plant from the chemistry to the plant structure. These calculations are grouped by scale and calculation type in different "methods".

##Inputs

###Required

1. Chemical Process Safety Sheet

  What: Basic safety properties of all the chemicals that appear in the design.
  
  Where: SDS sheets, Literature, Online databases
  
  Format: .txt
  
2. Design information report

  What: Process unit information for process conditions and functions
  
  Where: A simulator (e.g. Aspen Plus) or Rough hand calculations
  
  Format: .txt / .xls
  
3. Stream Table

  What: A table of the properties of all the streams in the process.
  
  Where: A simulator (e.g. Aspen Plus) / Rough hand calculations
  
  Format: .xls
  
4. Chemical Reactivity Worksheet (CRW) Report

  What: A government program that simulates a solution to tell if there will be any evolved gases or side reactions
  
  Where: The report of the program with all the chemicals in one solution
  
  Format: .txt
  

###Additional documents

1. Unit Alternatives

  What: A list of differences between the actual design and what is inputed in the design information report. There could be anything from end of the pipe technologies, to process control, to physical mitigation, to novel equipment.
  
  Where: A filled out rubric looking at the addition using a semi-quanitative risk ranking with three categories: Strength, Maintainability, and Reliability.
  
  Format: .xls
  
2. QSAR prediction reports (TEST)

  What: Predictions of chemical properties from their molecular structure. Many techniques can be used but for this procedure, USEPA's software TEST is automatically parsed.
  
  Where: The TEST output of a batch run of all the chemicals in the plant.
  
  Format: .txt

###IEHS specific databases
1. Failure rate per ASPEN block
2. Classification of ASPEN blocks


##Code notes


There is a config file with the location of all the inputs and the type of analysis that is loaded in "grabConfig". The "conceptual" design is the input-output structure. There is a alternative design out put in a spreadsheet for required values. The stream is the same format as well as the chemicals.

Most of the heavy assessment is done on "calcIntVal" which does all the assessment for the a single iteration including "Chem", "Mix", and "Vessel." Then the program iterates on the design with different inputs. In "calcHotSpot" it changes the conditions assessed over the length of the process unit. In "calcUUvals" , different conditions are fed into "calcIntVal" to simulate upsets. At the end of the analysis, "makeXlxSheets" summarizes the output into sheets with all the data. It is also available in the data structure like this: <br>

u(unitnumber, conditionNumber).method{agroType}(parameter). <br>

In other words, its in a array of structures with fields of the method names. The field itself is a cell of vectors as long as the aggregation possibilities. The final vector is as long as the number of parameters which is 13 for now.
