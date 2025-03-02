# Notes & Review of Report

### General
- your graphs look way nicer
- more organized code

### A1.2
specifically asks for a description in the task, just add my description, if you want:

"""
Apart from a huge jump between 2021 and 2022, the data follows vaguely a linear trendline, especially between 2018 and beginning 2021.
After the mentioned jump, the linear trend continues for approx. 6 motnhs in 2021. Then the model shows signs of stagnating growth for the total number of cars.
From end 2021 to 2024, the data follows almost a linear trend again, but with a smaller slope. 
This can be interpreted as a saturation in car registration across Denmark. Reasons could be: 
1. most of the population already has a car and no further need to register a new car
2. due to some legislation or other factor it has become less popular of viable to register a new car, which causes weaker demand for new car registrations
"""

### A2
- don't say 'sloppy' notation, its perfectly fine

### A2.1&2.2
- really good explanation and detail!!!

### A2.3&2.4
- we do get different values for the prediction intervals
- the graphs look similar though

### A3
first and fourth paragraph:
- the notion of "locality" might be a bit blurred by the other explanations
- it is that you basically null out certain values of the time-series. That is, how I understood, it's called "local" because you only consider a fraction of the global series
- write something about it in *3.1*, which you can C&P if you want

otherwise, very nice!

### A3.2
- maybe add in a weighted scatter plot of the data points :) I found that really illustrative.
- you can copy it from my notebook

### A3.3
- maybe note that you use the geometric series formula

### A3.4&5
- we get different parameter values (although I am pretty sure we already compared them and agreed), since we also get the same regression lines
- update: we do have the same values, but only in the notebook, not the PDF document yet
- you just have much bigger confidence intervals:
    - I think in the code you split up the equation that she provided in the R script for the intervals
    - you multiplied the sigma2 into the covariance matrix, I think it should be multiplied into the projection:

    ```
    covariance = sigma2 * np.linalg.inv(X.T @ W @ X)
    covariance_forecast = sigma2 + X_forecast @ covariance @ X_forecast.T

    #i think it should have been:
    covariance = np.linalg.inv(X.T @ W @ X)
    covariance_forecast = sigma2 + sigma * X_forecast @ covariance @ X_forecast.T
    ```


### A4.1
I uploaded the PDF into the folder

### A4.2
I actually did some more explanation for the "intuitive" equations, maybe you can copy in some of it

### A4.3
Figure 15: 
    1. where is the OLS curve?
    2. Is that a point-wise/step-wise RLS estimate?

### A4.4
Maybe copy in my multi-line plot with the vanishing regression lines; I think its actually a good plot to illustrate burn in and the effect of recursive forgetting

### A4.6
The plot with dots as min/max markers