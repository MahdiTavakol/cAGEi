#Developed by Mahdi Tavakol mahditavakol90@gmail.com
cAGEi:
	g++ main.cpp Atoms.cpp XST.cpp -o cAGEi
main:
	g++ -c main.cpp
Atoms:
	g++ -c Atoms.cpp
XST:
	g++ -c XST.cpp
clean:
	rm -rf *o cAGEi
	
