from preprocessing.run_preprocessing import run_preprocessing 
from hierarchical_aggregation.run_aggregation import run_aggregation 

def main():
    print("=== Starting full pipeline ===")

    run_preprocessing()
    run_aggregation()

    print("=== Pipeline complete ===")

if __name__ == "__main__":
    main()