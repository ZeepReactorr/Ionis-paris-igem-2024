# Aphidisperse - User Guide and Getting Started

## The problem

The main vector for the Beet Yellows Virus (BYV) is the aphid Myzus persicae. Cap'siRNA offers a targetted treatment to this virus, with maximized efficiency during the early stages of plant contamination. The best case scenario would therefore be to treat the beets a few days prior aphids invading the field, preventing them from transmitting the virus initially. The problematic here is therefore to be able to predict the time of arrival of the aphids in beet plantations, which is what we attempt to do here.

## Solution principle

We combine population dynamics, agent-based model and lattice-visualization to calculate the diffusion of aphids throughout the year in the French territory in order to have a rough appreciation of its migration behaviour using lattitude and longitude as its input. Those coordinates could very well be the coordinates of a farmer's culture for example, giving him a rough idea of when to apply pesticides. Aphidisperse was developped <b>from scratch</b> and does not rely on any external program to run its simulation. It is command line based, however a GUI is currently being developped to provide an easier use.

## Prerequisite

You will need to install ffmpeg in the Aphidisperse folder for the program to create the video. You can download the corresponding release here : https://ffmpeg.org/download.html

## Getting started

First clone the repository
```
git clone https://gitlab.igem.org/2024/software-tools/ionis-paris/
```

Move into the SafeRNA directory
```
cd ionis-paris/Aphidisperse
```

Install the necessary packages
```
pip install -r requirements.txt
```

## Use the program

```
python ~/PATH/TO/main.py lat lon year
```

Replace `lat` and `lon` by the coordinates you want. Currently coordinates oustide of France are not supported. The output of the program can be found in `~/PATH/TO/CLONING/DIRECTORY/output`. In the output you can find the starting state of the simulation for the desired year, as well as the video of the migration modeled. 


