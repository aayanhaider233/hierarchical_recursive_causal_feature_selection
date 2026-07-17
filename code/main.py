from preprocessing.run_preprocessing import run_preprocessing 

def main():
    print("=== Starting full pipeline ===")

    run_preprocessing()

    print("=== Pipeline complete ===")

if __name__ == "__main__":
    main()