import pandas as pd
import os
import matplotlib.pyplot as plt

def plot_scatter_from_csv(file_path, column_name):
    """
    Plot a scatter plot from a CSV file.

    Args:
        file_path (str): Path to the CSV file.
        column_name (str): Name of the column to be plotted.
    """
    # Read the CSV file into a pandas DataFrame
    predictions_df = pd.read_csv(file_path)

    # Generate x values (indices)
    x = range(len(predictions_df[column_name]))

    # Plotting the scatter plot
    plt.scatter(x, predictions_df[column_name], color='skyblue', edgecolor='black')

    # Adding labels and title
    plt.xlabel('Index')
    plt.ylabel('Predicted Values')
    plt.title('Scatter Plot of Predicted Values')

    # Showing plot
    plt.show()

