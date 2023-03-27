# Write Up Unscented Kalman Filter Project

The code mostly follows the example given in the lecture assignments. We follow the basic outline of the unscented Kalman Filter:

1. Initialize variables and constant matrices (like radar and lidar covariance) in the constructor. These probably could be marked as const or even constexpr and be initialized directly in the class definition, but were kept non-const to keep the initialization in the constructor
    1. Initial covariance for velocity and yaw angle are set high to avoid some possible instabilities
2. On first call of ProcessMeasurement the state vector x and the time stamp are initialized
3. On each subsequent call the object's state is predicted
    1. First the augmented sigma points are computed from the augmented state vector and augmented state covariance
    2. New sigma points are predicted according to the process model
    3. Using the sigma points the new state vector and covariance matrix can be calculated
4. Using the measurement package the belief on the object's state is updated
   1. The sigma points are transformed into the corresponding measurement space (Lidar or Radar) using the respective measurement model
   2. The innovation covariance matrix S is calculated and the independent noise matrix R_radar or R_lidar is added
   3. Kalman gain is computed from the cross correlation matrix (between state differences and residuals), weights for the sigma points are taken into account
   4. Object state is updated
   5. NIS is calculated to gauge quality of the chosen parameters

## Choice of Parameters

The Parameters were chosen such that the 5% NIS threshold (7.815 for radar, 5.991 for lidar) is exceeded about 5-10% of a time. I think this a decent indicator that the values for the process noise are chosen realistically, although a deeper statistical analysis of NIS could allow for a better parameter tuning.

## Short note on style

I tried to remove all explicit references to indices by moving the indices into enums. This makes the code safer and hopefully more readable. The main reason to do this, is to keep check of the correct State-Index correspondence. There were two minor issues with this approach: 1. I used C style enums instead of enum class to avoid the annoying cast to int when used as an index, 2. In older version of Eigen enum values cannot be used to construct vectors and matrices (although the underlying type is int), so I kept the state dimension as separate variables in the class and treated them as const.
