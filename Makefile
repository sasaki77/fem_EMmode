WORK_DIR = $(shell pwd)
BIN_DIR  = $(WORK_DIR)/bin


all:
	@if [ ! -d $(BIN_DIR) ]; \
		then echo "mkdir -p $(BIN_DIR)"; mkdir -p $(BIN_DIR); \
		fi
	cd $(WORK_DIR)/fem; make; cp fem $(BIN_DIR)
	cd $(WORK_DIR)/femvis; make; cp femvis $(BIN_DIR)

clean:
	cd $(WORK_DIR)/fem; make clean;
	cd $(WORK_DIR)/femvis; make clean;
