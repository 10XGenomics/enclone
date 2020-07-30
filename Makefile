# Build Enclone.
all: build

.DELETE_ON_ERROR:
.PHONY: all build test

# Build Enclone.
build:
	cargo build

# Test Enclone.
test: enclone_main/test/.gitignore
	git -C $(<D) fetch --depth=1 origin tag v0.4.45
	git -C $(<D) switch --detach FETCH_HEAD
	cargo test

# Download the test data.
enclone_main/test/.gitignore:
	git clone --depth=1 https://github.com/10XGenomics/enclone-data $(@D)
