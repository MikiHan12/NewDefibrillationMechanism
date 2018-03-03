# NewDefibrillationMechanism
Mechanism of Vortices Termination in Cardiac Muscle

## Description
Develop new, low-energy defribrillation methods that could be less damaging and less traumatic for the patient, and would save battery energy. However, these methods have not been entirely successful, due in part to an incomplete understanding of all the mechanisms present that may help or hinder the process of terminating the rotating waves present during fibrillation. Here we describe new mechanisms whereby a far-field electric field pulse terminates unpinned waves that are rotating in the vicinity of a blood vessel, plaque deposit or other heterogeneity in the gap junction conductivity.

## Method

Dimensionless Barkley equations is used to model the electrical activity of the heart. 
We run a series of two-dimentional computer simulations of a spiral wave rotating in the vicinity of a non-conducting obstacle. Then apply low-engergy electric field pulse for a short period of time. 
![picture](https://user-images.githubusercontent.com/25389100/36929586-940b4364-1e47-11e8-9d61-77af73db0a7c.png)


## Instruction

Compile all .cpp files in Terminal:

$ g++ -O *.cpp -o run

$./run

Run ContourUPlot.m in Matlab to see the simulation 

The following simulation mimics what happens electrically during fibrillation of the heart. 
