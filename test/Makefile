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

# 清理生成的文件
clean:
	rm -f $(TARGET)

# 伪目标
.PHONY: all clean