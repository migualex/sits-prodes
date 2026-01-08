# Satellite Image Time Series for Prodes Amazônia Analysis

The sits package uses satellite image time series for land classification through a *time-first, space-later approach*.  
In the data preparation process, collections of large Earth observation images are organized into data cubes.  
Each spatial location in a data cube is associated with a time series. Locations with known labels are used to train a machine learning algorithm, which then classifies all time series within the cube, as illustrated below.

<p align="center">
  <img src="https://github.com/user-attachments/assets/65decf6f-3782-49c5-b608-a71d73c57ad0" alt="sits_general_view" width="80%">
</p>

<p align="center">
  <em>Figure 1: General view of the sits workflow.</em>
</p>

---

The sits API provides a set of functions that can be linked together to create an end-to-end workflow for land classification.  
At its core, the sits package offers eight main functions, as illustrated in the figure below:

1. **sits_cube()** – Extract data from an analysis-ready data (ARD) collection, producing a *non-regular data cube* object.  
2. **sits_regularize()** – Convert a non-regular data cube into a *regular* one, required for training machine learning algorithms.  
3. **sits_apply()** – Obtain new bands and indices by applying operations to regular data cubes.  
4. **sits_get_data()** – Using ground truth values (e.g., CSV or SHP) and a regular data cube, extract training samples containing time series for selected locations.  
5. **sits_train()** – Select a machine learning algorithm and train a classification model.  
6. **sits_classify()** – Apply the trained model to classify a regular data cube, generating a *probability cube* with class probabilities for each pixel.  
7. **sits_smooth()** – Remove outliers from the probability cube.  
8. **sits_label_classification()** – Generate a *thematic map* from a smoothed probability cube.

<p align="center">
  <img src="https://github.com/user-attachments/assets/b1d616f4-8273-4046-b1b9-b8a09ef15d22" alt="sits_api" width="85%">
</p>

<p align="center">
  <em>Figure 2: Main functions of the sits API.</em>
</p>

---

Special thanks to the Brazil Data Cube team for their support and training.

_Source: Rolf Simoes, Gilberto Camara, Gilberto Queiroz, Felipe Souza, Pedro R. Andrade, Lorena Santos, Alexandre Carvalho, and Karine Ferreira. Satellite Image Time Series Analysis for Big Earth Observation Data. Remote Sensing, 13, p. 2428, 2021._
