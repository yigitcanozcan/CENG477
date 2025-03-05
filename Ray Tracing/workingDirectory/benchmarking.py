import os
from timeit import default_timer
import json

inputList = os.listdir("inputs")
times = dict()



for input in inputList:
    start = default_timer()
    os.system("./raytracer {}".format(os.path.join("inputs",input)))
    end = default_timer()
    duration = end-start
    times[input] = duration
    print("{} took {} seconds to render".format(input,duration))

with open("benchmark.json","w") as f:
    json.dump(times,f)   
