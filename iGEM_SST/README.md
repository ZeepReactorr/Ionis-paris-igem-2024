# iGEM team Subject Search Tool

## The problem

In the early stages of our brainstorming for our project, a recurrent questions was "but what did the other teams do before?". And when we were finding a good idea, we would ask "Has it already been done in iGEM?" or "Are there project on which we can base ourselves on?". And we could not find efficiently those teams, because there was no keyword search tool available to browse easily the previous iGEM teams.

## The solution 

Well, then we decided to make our own tool! First, we scrapped the freely available data from all the teams since the beginning of iGEM available in the team list of the official website https://teams.igem.org/. The details collected were : 
- Wiki url
- Village
- Abstract if available

We then developped a keyword-based program, able to take in several keywords, and returning in a file (if you use the command line version of the tool), or directly on your screen (if you use our streamlit webapp or our stand-alone application), all the teams with the desired keyword in their abstract!

## Getting started

If you want to use the command line version of the tool, clone the repository :
```
git clone https://gitlab.igem.org/2024/software-tools/ionis-paris/
```

This project is based on the stdlib of Python, so no additional package are required. To use the tool from the command line interface just type :
```
python ~/PATH/TO/GUI.py
```

A rudimentary GUI will appear, with a search bar. Just type your keywords, and voilà! the teams will be displayed!

Another easy way is to access the tool through the streamlit release, available at this link : https://igemteamsearch.streamlit.app/ <br>
The app may go down due to not being used often (for now!).