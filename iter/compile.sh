g++ -std=c++17 -O2 banemecp.cpp -o banemecp -lstdc++fs
g++ -std=c++17 -O3 -DNDEBUG -static-libgcc -static-libstdc++ -s banemecp.cpp -o banemecp -lstdc++fs