# 定义编译器
CXX = g++

# 定义源文件和目标文件
SRC = test.cpp
TARGET = test

# 默认目标
all: $(TARGET)

# 链接源文件生成可执行文件
$(TARGET): $(SRC)
	$(CXX) -o $@ $^

draw: draw.cpp
	$(CXX) draw.cpp -I/usr/include/python3.10 -I~/.local/lib/python3.10/site-packages/numpy/core/include -lpython3.10 -o draw
curve: curve.cpp
	$(CXX) curve.cpp -I/usr/include/python3.10 -I~/.local/lib/python3.10/site-packages/numpy/core/include -lpython3.10 -o curve
cycle: cycle.cpp
	$(CXX) cycle.cpp -I/usr/include/python3.10 -I~/.local/lib/python3.10/site-packages/numpy/core/include -lpython3.10 -o cycle
# 清理生成的文件
clean:
	rm -f draw
	rm -f $(TARGET)
	rm -f curve
	rm -f testBC
	rm -f testmap
	rm -f testSP
	rm -f testSplines
	rm -f testVec
	rm -f testBisection
	rm -f cycle

# 伪目标
.PHONY: all clean