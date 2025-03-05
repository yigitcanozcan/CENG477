g++ rayTracer.cpp -o rayTracer
cat sceneDescription.txt | ./rayTracer > output.ppm
display output.ppm
echo 'DONE'
