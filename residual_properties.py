import pandas as pd

# Data for the H^R0
data_HR0 = [
    [-6.045, -6.043, -6.040, -6.034, -6.022, -6.011, -5.999, -5.987, -5.975, -5.957, -5.927, -5.868, -5.748, -5.628, -5.446],
    [-5.906, -5.904, -5.901, -5.895, -5.882, -5.870, -5.858, -5.845, -5.833, -5.814, -5.783, -5.721, -5.595, -5.469, -5.278],
    [-5.763, -5.761, -5.757, -5.751, -5.738, -5.726, -5.713, -5.700, -5.687, -5.668, -5.636, -5.572, -5.442, -5.311, -5.113],
    [-5.615, -5.612, -5.609, -5.603, -5.590, -5.577, -5.564, -5.551, -5.538, -5.519, -5.486, -5.421, -5.288, -5.154, -5.950],
    [-5.465, -5.463, -5.459, -5.453, -5.440, -5.427, -5.414, -5.401, -5.388, -5.369, -5.336, -5.279, -5.135, -4.999, -4.791],
    [-0.032, -5.312, -5.309, -5.303, -5.290, -5.278, -5.265, -5.252, -5.239, -5.220, -5.187, -5.121, -4.986, -4.849, -4.638],
    [-0.027, -5.162, -5.159, -5.153, -5.141, -5.129, -5.116, -5.104, -5.091, -5.073, -5.041, -4.976, -4.842, -4.794, -4.492],
    [-0.023, -0.118, -5.008, -5.002, -4.991, -4.980, -4.968, -4.956, -4.949, -4.927, -4.896, -4.833, -4.702, -4.565, -4.353],
    [-0.020, -0.101, -0.213, -4.848, -4.838, -4.828, -4.818, -4.808, -4.797, -4.781, -4.752, -4.693, -4.566, -4.432, -4.221],
    [-0.017, -0.088, -0.183, -4.687, -4.679, -4.672, -4.664, -4.655, -4.646, -4.632, -4.607, -4.554, -4.434, -4.393, -4.095],
    [-0.015, -0.078, -0.160, -0.345, -4.507, -4.504, -4.499, -4.494, -4.488, -4.478, -4.459, -4.413, -4.303, -4.178, -3.974],
    [-0.014, -0.069, -0.141, -0.300, -4.309, -4.313, -4.316, -4.316, -4.316, -4.312, -4.302, -4.269, -4.173, -4.056, -3.857],
    [-0.012, -0.062, -0.126, -0.264, -0.596, -4.074, -4.094, -4.108, -4.118, -4.127, -4.132, -4.119, -4.043, -3.935, -3.744],
    [-0.011, -0.058, -0.118, -0.246, -0.545, -0.960, -3.920, -3.953, -3.976, -4.000, -4.020, -4.024, -3.963, -3.863, -3.678],
    [-0.011, -0.056, -0.113, -0.235, -0.516, -0.885, -3.763, -3.825, -3.865, -3.904, -3.940, -3.958, -3.910, -3.815, -3.634],
    [-0.011, -0.054, -0.109, -0.225, -0.490, -0.824, -1.356, -3.658, -3.732, -3.796, -3.853, -3.890, -3.856, -3.767, -3.591],
    [-0.010, -0.053, -0.107, -0.221, -0.478, -0.797, -1.273, -3.544, -3.652, -3.736, -3.806, -3.854, -3.829, -3.743, -3.569],
    [-0.010, -0.052, -0.105, -0.216, -0.466, -0.773, -1.206, -3.376, -3.558, -3.670, -3.758, -3.818, -3.801, -3.719, -3.548],
    [-0.010, -0.051, -0.103, -0.212, -0.455, -0.750, -1.151, -2.584, -3.441, -3.598, -3.706, -3.782, -3.774, -3.695, -3.526],
    [-0.010, -0.050, -0.101, -0.208, -0.445, -0.721, -1.102, -1.796, -3.283, -3.516, -3.652, -3.744, -3.746, -3.671, -3.505],
    [-0.010, -0.049, -0.099, -0.203, -0.434, -0.708, -1.060, -1.627, -3.039, -3.422, -3.595, -3.705, -3.718, -3.647, -3.484],
    [-0.009, -0.046, -0.094, -0.192, -0.407, -0.654, -0.955, -1.359, -2.034, -3.030, -3.398, -3.583, -3.632, -3.575, -3.420],
    [-0.008, -0.042, -0.086, -0.175, -0.367, -0.581, -0.827, -1.120, -1.487, -2.203, -2.965, -3.353, -3.484, -3.453, -3.315],
    [-0.008, -0.039, -0.079, -0.160, -0.334, -0.523, -0.732, -0.968, -1.239, -1.719, -2.479, -3.091, -3.329, -3.329, -3.211],
    [-0.007, -0.036, -0.073, -0.148, -0.305, -0.474, -0.657, -0.857, -1.076, -1.443, -2.079, -2.801, -3.166, -3.202, -3.107],
    [-0.006, -0.031, -0.063, -0.127, -0.259, -0.399, -0.545, -0.698, -0.860, -1.116, -1.560, -2.274, -2.825, -2.942, -2.899],
    [-0.005, -0.027, -0.055, -0.110, -0.224, -0.341, -0.463, -0.588, -0.716, -0.915, -1.253, -1.857, -2.486, -2.679, -2.692],
    [-0.005, -0.024, -0.048, -0.097, -0.196, -0.297, -0.400, -0.505, -0.611, -0.774, -1.046, -1.549, -2.175, -2.421, -2.486],
    [-0.004, -0.021, -0.043, -0.086, -0.173, -0.261, -0.350, -0.440, -0.531, -0.667, -0.894, -1.318, -1.904, -2.177, -2.285],
    [-0.004, -0.019, -0.038, -0.076, -0.153, -0.231, -0.309, -0.387, -0.446, -0.583, -0.777, -1.139, -1.672, -1.953, -2.091],
    [-0.003, -0.017, -0.034, -0.068, -0.137, -0.206, -0.275, -0.344, -0.413, -0.515, -0.683, -0.996, -1.476, -1.751, -1.908],
    [-0.003, -0.015, -0.031, -0.062, -0.123, -0.185, -0.246, -0.307, -0.368, -0.458, -0.606, -0.880, -1.309, -1.571, -1.736],
    [-0.003, -0.014, -0.028, -0.056, -0.111, -0.167, -0.222, -0.276, -0.330, -0.411, -0.541, -0.782, -1.167, -1.411, -1.577],
    [-0.002, -0.012, -0.023, -0.046, -0.092, -0.137, -0.182, -0.226, -0.269, -0.334, -0.437, -0.629, -0.937, -1.143, -1.295],
    [-0.002, -0.010, -0.019, -0.038, -0.076, -0.114, -0.150, -0.187, -0.222, -0.275, -0.359, -0.513, -0.761, -0.929, -1.058],
    [-0.002, -0.008, -0.016, -0.032, -0.064, -0.095, -0.125, -0.155, -0.185, -0.228, -0.297, -0.422, -0.621, -0.756, -0.858],
    [-0.001, -0.007, -0.014, -0.027, -0.054, -0.080, -0.105, -0.130, -0.154, -0.190, -0.246, -0.348, -0.508, -0.614, -0.689],
    [-0.001, -0.006, -0.011, -0.023, -0.045, -0.067, -0.088, -0.109, -0.129, -0.159, -0.205, -0.288, -0.415, -0.495, -0.545],
    [-0.001, -0.004, -0.007, -0.015, -0.029, -0.043, -0.056, -0.069, -0.081, -0.099, -0.127, -0.174, -0.239, -0.270, -0.264],
    [-0.000, -0.002, -0.005, -0.009, -0.017, -0.026, -0.033, -0.041, -0.048, -0.058, -0.072, -0.095, -0.116, -0.110, -0.061]
]

# Headers and index
header = [0.0100, 0.0500, 0.1000, 0.2000, 0.4000, 0.6000, 0.8000, 1.0000, 1.2000, 1.5000, 2.0000, 3.0000, 5.0000, 7.0000, 10.000]
index = [0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.93, 0.95, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00]

HR0 = pd.DataFrame(data_HR0, columns=header, index=index)


# Change the name of the header and index
HR0 = HR0.rename_axis("Pr", axis=1)
HR0 = HR0.rename_axis("Tr", axis=0)


# Data for the H^R1
data_HR1 = [
    [-11.098, -11.096, -11.095, -11.091, -11.083, -11.076, -11.069, -11.062, -11.055, -11.044, -11.027, -10.992, -10.935, -10.872, -10.781],
    [-10.656, -10.655, -10.654, -10.653, -10.650, -10.646, -10.643, -10.640, -10.637, -10.632, -10.624, -10.609, -10.581, -10.554, -10.529],
    [-10.121, -10.121, -10.121, -10.120, -10.121, -10.121, -10.121, -10.121, -10.121, -10.121, -10.122, -10.123, -10.128, -10.135, -10.150],
    [-9.515, -9.515, -9.516, -9.517, -9.519, -9.521, -9.523, -9.525, -9.527, -9.531, -9.537, -9.549, -9.576, -9.611, -9.663],
    [-8.868, -8.869, -8.870, -8.872, -8.876, -8.880, -8.884, -8.888, -8.892, -8.899, -8.909, -8.932, -8.978, -9.030, -9.111],
    [-0.080, -8.211, -8.212, -8.215, -8.221, -8.226, -8.232, -8.238, -8.243, -8.243, -8.252, -8.267, -8.298, -8.360, -8.245],
    [-0.059, -7.568, -7.570, -7.573, -7.579, -7.585, -7.591, -7.596, -7.603, -7.611, -7.632, -7.669, -7.745, -7.824, -7.950],
    [-0.045, -0.247, -6.949, -6.952, -6.959, -6.966, -6.973, -6.980, -6.987, -6.997, -7.017, -7.059, -7.147, -7.239, -7.381],
    [-0.034, -0.185, -0.415, -6.360, -6.367, -6.373, -6.381, -6.388, -6.395, -6.407, -6.429, -6.475, -6.574, -6.677, -6.837],
    [-0.027, -0.142, -0.306, -5.796, -5.802, -5.809, -5.816, -5.824, -5.832, -5.845, -5.868, -5.918, -6.027, -6.142, -6.318],
    [-0.021, -0.110, -0.234, -0.542, -5.266, -5.271, -5.278, -5.285, -5.293, -5.306, -5.330, -5.385, -5.506, -5.632, -5.824],
    [-0.017, -0.087, -0.182, -0.401, -4.753, -4.754, -4.758, -4.763, -4.771, -4.784, -4.810, -4.872, -5.000, -5.149, -5.358],
    [-0.014, -0.070, -0.144, -0.308, -0.751, -4.254, -4.248, -4.249, -4.255, -4.268, -4.298, -4.371, -4.530, -4.688, -4.916],
    [-0.012, -0.061, -0.126, -0.265, -0.612, -1.236, -3.942, -3.934, -3.937, -3.951, -3.987, -4.073, -4.251, -4.422, -4.662],
    [-0.011, -0.056, -0.115, -0.241, -0.542, -0.994, -3.737, -3.712, -3.713, -3.730, -3.773, -3.873, -4.068, -4.248, -4.497],
    [-0.010, -0.052, -0.105, -0.219, -0.483, -0.837, -1.616, -3.470, -3.467, -3.492, -3.551, -3.670, -3.885, -4.077, -4.336],
    [-0.010, -0.050, -0.101, -0.209, -0.457, -0.776, -1.324, -3.332, -3.327, -3.363, -3.434, -3.568, -3.795, -3.992, -4.257],
    [-0.009, -0.048, -0.097, -0.200, -0.433, -0.722, -1.154, -3.164, -3.164, -3.223, -3.313, -3.464, -3.705, -3.909, -4.178],
    [-0.009, -0.046, -0.093, -0.191, -0.410, -0.675, -1.034, -2.471, -2.952, -3.065, -3.186, -3.358, -3.615, -3.825, -4.100],
    [-0.009, -0.044, -0.089, -0.183, -0.389, -0.632, -0.940, -1.375, -2.595, -2.880, -3.051, -3.251, -3.525, -3.742, -4.023],
    [-0.008, -0.042, -0.085, -0.175, -0.370, -0.594, -0.863, -1.180, -1.723, -2.650, -2.906, -3.142, -3.435, -3.661, -3.947],
    [-0.007, -0.037, -0.075, -0.153, -0.318, -0.498, -0.691, -0.877, -0.878, -1.496, -2.381, -2.800, -3.167, -3.418, -3.722],
    [-0.006, -0.030, -0.061, -0.123, -0.251, -0.381, -0.507, -0.617, -0.673, -0.617, -1.261, -2.167, -2.720, -3.023, -3.362],
    [-0.005, -0.025, -0.050, -0.099, -0.199, -0.296, -0.385, -0.459, -0.503, -0.487, -0.604, -1.497, -2.275, -2.641, -3.019],
    [-0.004, -0.020, -0.040, -0.080, -0.158, -0.232, -0.297, -0.349, -0.381, -0.381, -0.361, -0.934, -1.840, -2.273, -2.692],
    [-0.003, -0.013, -0.026, -0.052, -0.100, -0.142, -0.177, -0.203, -0.218, -0.218, -0.178, -0.300, -1.066, -1.592, -2.086],
    [-0.002, -0.008, -0.016, -0.032, -0.060, -0.083, -0.100, -0.111, -0.115, -0.128, -0.070, -0.044, -0.504, -1.012, -1.547],
    [-0.001, -0.005, -0.009, -0.018, -0.032, -0.042, -0.048, -0.049, -0.046, -0.032, 0.008, 0.078, -0.142, -0.556, -1.080],
    [0.000, -0.002, -0.004, -0.007, -0.012, -0.013, -0.011, -0.005, 0.004, 0.023, 0.065, 0.151, 0.082, -0.217, -0.689],
    [0.000, 0.000, 0.000, 0.000, 0.003, 0.009, 0.017, 0.027, 0.040, 0.063, 0.109, 0.202, 0.223, 0.028, -0.369],
    [0.000, 0.001, 0.003, 0.006, 0.015, 0.025, 0.037, 0.051, 0.067, 0.094, 0.143, 0.241, 0.317, 0.203, -0.112],
    [0.001, 0.003, 0.005, 0.011, 0.023, 0.037, 0.053, 0.070, 0.088, 0.117, 0.169, 0.271, 0.381, 0.330, 0.092],
    [0.001, 0.003, 0.007, 0.015, 0.030, 0.047, 0.065, 0.085, 0.105, 0.136, 0.190, 0.295, 0.428, 0.424, 0.255],
    [0.001, 0.005, 0.010, 0.020, 0.040, 0.062, 0.083, 0.106, 0.128, 0.163, 0.221, 0.331, 0.493, 0.551, 0.489],
    [0.001, 0.006, 0.012, 0.023, 0.047, 0.071, 0.095, 0.120, 0.144, 0.181, 0.242, 0.356, 0.535, 0.631, 0.645],
    [0.001, 0.006, 0.013, 0.026, 0.052, 0.078, 0.104, 0.130, 0.156, 0.194, 0.257, 0.376, 0.567, 0.687, 0.754],
    [0.001, 0.007, 0.014, 0.028, 0.055, 0.082, 0.110, 0.137, 0.164, 0.204, 0.269, 0.391, 0.591, 0.729, 0.836],
    [0.001, 0.007, 0.014, 0.029, 0.058, 0.086, 0.114, 0.142, 0.170, 0.211, 0.278, 0.403, 0.611, 0.763, 0.899],
    [0.002, 0.008, 0.016, 0.031, 0.062, 0.092, 0.122, 0.152, 0.181, 0.224, 0.294, 0.425, 0.650, 0.827, 1.015],
    [0.002, 0.008, 0.016, 0.032, 0.064, 0.096, 0.127, 0.158, 0.188, 0.233, 0.306, 0.442, 0.680, 0.874, 1.097]
]


HR1 = pd.DataFrame(data_HR1, columns=header, index=index)


# Change the name of the header and index
HR1 = HR1.rename_axis("Pr", axis=1)
HR1 = HR1.rename_axis("Tr", axis=0)


data_SR0 = [
    [-11.614, -10.008, -9.319, -8.635, -7.961, -7.574, -7.304, -7.099, -6.935, -6.740, -6.497, -6.180, -5.847, -5.683, -5.578],
    [-11.185, -9.579, -8.890, -8.205, -7.529, -7.140, -6.869, -6.663, -6.497, -6.299, -6.052, -5.728, -5.376, -5.194, -5.060],
    [-10.802, -9.196, -8.506, -7.821, -7.144, -6.755, -6.483, -6.275, -6.109, -5.909, -5.660, -5.330, -4.967, -4.772, -4.619],
    [-10.453, -8.847, -8.157, -7.472, -6.794, -6.404, -6.132, -5.924, -5.757, -5.557, -5.306, -4.974, -4.603, -4.401, -4.234],
    [-10.137, -8.531, -7.841, -7.156, -6.479, -6.089, -5.816, -5.608, -5.441, -5.240, -4.989, -4.656, -4.282, -4.074, -3.899],
    [-0.038, -8.245, -7.555, -6.870, -6.193, -5.803, -5.531, -5.324, -5.157, -4.956, -4.706, -4.373, -3.998, -3.788, -3.607],
    [-0.029, -7.983, -7.294, -6.610, -5.933, -5.544, -5.273, -5.066, -4.900, -4.700, -4.451, -4.120, -3.747, -3.537, -3.353],
    [-0.023, -0.122, -7.052, -6.368, -5.694, -5.306, -5.036, -4.830, -4.665, -4.467, -4.220, -3.892, -3.523, -3.315, -3.131],
    [-0.018, -0.096, -0.206, -6.140, -5.467, -5.082, -4.814, -4.610, -4.446, -4.250, -4.007, -3.684, -3.322, -3.117, -2.935],
    [-0.015, -0.078, -0.164, -5.917, -5.248, -4.866, -4.600, -4.399, -4.238, -4.045, -3.807, -3.491, -3.138, -2.939, -2.761],
    [-0.013, -0.064, -0.134, -0.294, -5.026, -4.694, -4.388, -4.191, -4.034, -3.846, -3.615, -3.310, -2.970, -2.777, -2.605],
    [-0.011, -0.054, -0.111, -0.239, -4.785, -4.418, -4.166, -3.976, -3.825, -3.646, -3.425, -3.135, -2.812, -2.629, -2.463],
    [-0.009, -0.046, -0.094, -0.199, -0.463, -4.145, -3.912, -3.738, -3.599, -3.434, -3.231, -2.964, -2.663, -2.491, -2.334],
    [-0.008, -0.042, -0.085, -0.179, -0.408, -0.750, -3.723, -3.569, -3.444, -3.295, -3.108, -2.860, -2.577, -2.412, -2.262],
    [-0.008, -0.039, -0.080, -0.168, -0.377, -0.671, -3.556, -3.433, -3.326, -3.193, -3.023, -2.790, -2.520, -2.362, -2.215],
    [-0.007, -0.037, -0.075, -0.157, -0.350, -0.607, -1.056, -3.259, -3.188, -3.081, -2.932, -2.719, -2.463, -2.312, -2.170],
    [-0.007, -0.036, -0.073, -0.153, -0.337, -0.580, -0.971, -3.142, -3.106, -3.019, -2.884, -2.682, -2.436, -2.287, -2.148],
    [-0.007, -0.035, -0.071, -0.148, -0.326, -0.555, -0.903, -2.972, -3.010, -2.953, -2.835, -2.646, -2.408, -2.263, -2.126],
    [-0.007, -0.034, -0.069, -0.144, -0.315, -0.532, -0.847, -2.178, -2.893, -2.879, -2.784, -2.609, -2.380, -2.239, -2.105],
    [-0.007, -0.033, -0.067, -0.139, -0.304, -0.510, -0.799, -1.391, -2.736, -2.798, -2.730, -2.571, -2.352, -2.215, -2.083],
    [-0.006, -0.032, -0.065, -0.135, -0.294, -0.491, -0.757, -1.225, -2.495, -2.706, -2.673, -2.533, -2.325, -2.191, -2.062],
    [-0.006, -0.030, -0.060, -0.124, -0.267, -0.439, -0.656, -0.965, -1.523, -2.328, -2.483, -2.415, -2.242, -2.121, -2.001],
    [-0.005, -0.026, -0.053, -0.108, -0.230, -0.371, -0.537, -0.742, -1.012, -1.557, -2.081, -2.202, -2.104, -2.007, -1.903],
    [-0.005, -0.023, -0.047, -0.096, -0.201, -0.319, -0.452, -0.607, -0.790, -1.126, -1.649, -1.968, -1.966, -1.897, -1.810],
    [-0.004, -0.021, -0.042, -0.085, -0.177, -0.277, -0.389, -0.512, -0.651, -0.890, -1.308, -1.727, -1.827, -1.789, -1.722],
    [-0.003, -0.017, -0.033, -0.068, -0.140, -0.217, -0.298, -0.385, -0.478, -0.628, -0.891, -1.299, -1.554, -1.581, -1.556],
    [-0.003, -0.014, -0.027, -0.056, -0.114, -0.174, -0.237, -0.303, -0.375, -0.478, -0.663, -0.990, -1.303, -1.386, -1.402],
    [-0.002, -0.011, -0.023, -0.046, -0.094, -0.143, -0.194, -0.246, -0.299, -0.381, -0.520, -0.777, -1.088, -1.208, -1.260],
    [-0.002, -0.010, -0.019, -0.039, -0.079, -0.120, -0.162, -0.204, -0.247, -0.312, -0.421, -0.628, -0.913, -1.050, -1.130],
    [-0.002, -0.008, -0.017, -0.033, -0.067, -0.102, -0.137, -0.172, -0.208, -0.261, -0.350, -0.519, -0.773, -0.915, -1.013],
    [-0.001, -0.007, -0.014, -0.029, -0.058, -0.088, -0.117, -0.147, -0.177, -0.222, -0.296, -0.438, -0.661, -0.799, -0.908],
    [-0.001, -0.006, -0.013, -0.025, -0.051, -0.076, -0.102, -0.127, -0.153, -0.191, -0.255, -0.375, -0.570, -0.702, -0.815],
    [-0.001, -0.006, -0.011, -0.022, -0.044, -0.067, -0.089, -0.111, -0.134, -0.167, -0.221, -0.625, -0.497, -0.620, -0.733],
    [-0.001, -0.004, -0.009, -0.018, -0.035, -0.053, -0.070, -0.087, -0.105, -0.130, -0.172, -0.251, -0.388, -0.492, -0.599],
    [-0.001, -0.004, -0.007, -0.014, -0.028, -0.042, -0.056, -0.070, -0.084, -0.104, -0.138, -0.201, -0.311, -0.399, -0.496],
    [-0.001, -0.003, -0.006, -0.012, -0.023, -0.035, -0.046, -0.058, -0.069, -0.086, -0.113, -0.164, -0.255, -0.329, -0.416],
    [0.000, -0.002, -0.005, -0.010, -0.020, -0.029, -0.039, -0.048, -0.058, -0.072, -0.094, -0.137, -0.213, -0.277, -0.353],
    [0.000, -0.002, -0.004, -0.008, -0.017, -0.025, -0.033, -0.041, -0.049, -0.061, -0.080, -0.116, -0.181, -0.236, -0.303],
    [0.000, -0.001, -0.003, -0.006, -0.012, -0.017, -0.023, -0.029, -0.034, -0.042, -0.056, -0.081, -0.126, -0.166, -0.216],
    [0.000, -0.001, -0.002, -0.004, -0.009, -0.013, -0.017, -0.021, -0.025, -0.031, -0.041, -0.059, -0.093, -0.123, -0.162]
]

SR0 = pd.DataFrame(data_SR0, columns=header, index=index)


# Change the name of the header and index
SR0 = SR0.rename_axis("Pr", axis=1)
SR0 = SR0.rename_axis("Tr", axis=0)


data_SR1 = [
    [ -16.782, -16.774, -16.764, -16.744, -16.705, -16.665, -16.626, -16.586, -16.547, -16.488, -16.390, -16.195, -15.837, -15.468, -14.925],
    [ -15.413, -15.408, -15.401, -15.387, -15.359, -15.333, -15.305, -15.278, -15.251, -15.211, -15.144, -15.011, -14.751, -14.496, -14.153],
    [ -13.990, -13.986, -13.981, -13.972, -13.953, -13.934, -13.915, -13.896, -13.877, -13.849, -13.803, -13.714, -13.541, -13.376, -13.144],
    [ -12.564, -12.561, -12.558, -12.551, -12.537, -12.523, -12.509, -12.496, -12.482, -12.462, -12.430, -12.367, -12.248, -12.145, -11.999],
    [ -11.202, -11.200, -11.197, -11.092, -11.082, -11.172, -11.162, -11.153, -11.143, -11.129, -11.107, -11.063, -10.985, -10.920, -10.836],
    [ -0.115, -9.948, -9.946, -9.942, -9.935, -9.928, -9.921, -9.914, -9.907, -9.897, -9.882, -9.853, -9.806, -9.769, -9.732],
    [ -0.078, -8.828, -8.826, -8.823, -8.817, -8.811, -8.806, -8.799, -8.794, -8.787, -8.777, -8.760, -8.736, -8.723, -8.720],
    [ -0.055, -0.309, -7.832, -7.829, -7.824, -7.819, -7.815, -7.510, -7.807, -7.801, -7.794, -7.784, -7.779, -7.785, -7.811],
    [ -0.040, -0.216, -0.491, -6.951, -6.945, -6.941, -6.937, -6.933, -6.930, -6.926, -6.922, -6.919, -6.929, -6.952, -7.002],
    [ -0.029, -0.156, -0.340, -6.173, -6.167, -6.162, -6.158, -6.155, -6.152, -6.149, -6.147, -6.149, -6.174, -6.213, -6.285],
    [ -0.022, -0.116, -0.246, -0.578, -5.475, -5.468, -5.462, -5.458, -5.455, -5.453, -5.452, -5.461, -5.501, -5.555, -5.648],
    [ -0.017, -0.088, -0.183, -0.400, -4.853, -4.841, -4.832, -4.826, -4.822, -4.820, -4.822, -4.839, -4.898, -4.969, -5.082],
    [ -0.013, -0.068, -0.140, -0.301, -0.744, -4.269, -4.249, -4.238, -4.232, -4.230, -4.236, -4.267, -4.351, -4.442, -4.578],
    [ -0.011, -0.058, -0.120, -0.254, -0.593, -1.219, -3.914, -3.894, -3.885, -3.884, -3.896, -3.941, -4.046, -4.151, -4.300],
    [ -0.010, -0.053, -0.109, -0.228, -0.517, -0.961, -3.697, -3.658, -3.647, -3.648, -3.669, -3.728, -3.851, -3.966, -4.125],
    [ -0.010, -0.048, -0.099, -0.206, -0.456, -0.797, -1.570, -3.406, -3.391, -3.401, -3.437, -3.517, -3.661, -3.788, -3.957],
    [ -0.009, -0.046, -0.094, -0.196, -0.429, -0.734, -1.270, -3.264, -3.247, -3.268, -3.318, -3.412, -3.569, -3.701, -3.875],
    [ -0.009, -0.044, -0.090, -0.186, -0.405, -0.680, -1.098, -3.093, -3.082, -3.126, -3.195, -3.306, -3.477, -3.616, -3.796],
    [ -0.008, -0.042, -0.086, -0.177, -0.382, -0.632, -0.977, -2.399, -2.868, -2.967, -3.067, -3.200, -3.387, -3.532, -3.717],
    [ -0.008, -0.040, -0.082, -0.169, -0.361, -0.590, -0.883, -1.306, -2.513, -2.784, -2.933, -3.094, -3.297, -3.450, -3.640],
    [ -0.008, -0.039, -0.078, -0.161, -0.342, -0.552, -0.807, -1.113, -1.655, -2.557, -2.790, -2.986, -3.209, -3.369, -3.565],
    [ -0.007, -0.034, -0.069, -0.140, -0.292, -0.460, -0.642, -0.820, -0.831, -1.443, -2.283, -2.655, -2.949, -3.134, -3.348],
    [ -0.005, -0.028, -0.055, -0.112, -0.229, -0.350, -0.470, -0.577, -0.640, -0.618, -1.241, -2.067, -2.534, -2.767, -3.013],
    [ -0.005, -0.023, -0.045, -0.091, -0.183, -0.275, -0.361, -0.437, -0.489, -0.502, -0.654, -1.471, -2.138, -2.428, -2.708],
    [ -0.004, -0.019, -0.037, -0.075, -0.149, -0.220, -0.286, -0.343, -0.385, -0.412, -0.447, -0.991, -1.767, -2.115, -2.430],
    [ -0.003, -0.013, -0.026, -0.052, -0.102, -0.148, -0.190, -0.226, -0.254, -0.282, -0.300, -0.481, -1.147, -1.569, -1.944],
    [ -0.002, -0.010, -0.019, -0.037, -0.072, -0.104, -0.133, -0.158, -0.178, -0.200, -0.220, -0.290, -0.730, -1.138, -1.544],
    [ -0.001, -0.007, -0.014, -0.027, -0.053, -0.076, -0.097, -0.115, -0.130, -0.147, -0.166, -0.206, -0.479, -0.823, -1.222],
    [ -0.001, -0.005, -0.011, -0.021, -0.040, -0.057, -0.073, -0.086, -0.098, -0.112, -0.129, -0.159, -0.334, -0.604, -0.969],
    [ -0.001, -0.004, -0.008, -0.016, -0.031, -0.044, -0.056, -0.067, -0.076, -0.087, -0.102, -0.127, -0.248, -0.456, -0.775],
    [ -0.001, -0.003, -0.006, -0.013, -0.024, -0.035, -0.044, -0.053, -0.060, -0.070, -0.083, -0.105, -0.195, -0.355, -0.628],
    [ -0.001, -0.003, -0.005, -0.010, -0.019, -0.028, -0.036, -0.043, -0.049, -0.057, -0.069, -0.089, -0.160, -0.286, -0.518],
    [ -0.000, -0.002, -0.004, -0.008, -0.016, -0.023, -0.029, -0.035, -0.040, -0.048, -0.058, -0.077, -0.136, -0.238, -0.434],
    [ -0.000, -0.001, -0.003, -0.006, -0.011, -0.016, -0.021, -0.025, -0.029, -0.035, -0.043, -0.060, -0.105, -0.178, -0.322],
    [ -0.000, -0.001, -0.002, -0.004, -0.008, -0.012, -0.015, -0.019, -0.022, -0.027, -0.034, -0.048, -0.086, -0.143, -0.254],
    [ -0.000, -0.001, -0.002, -0.003, -0.006, -0.009, -0.012, -0.015, -0.018, -0.021, -0.028, -0.041, -0.074, -0.120, -0.210],
    [ -0.000, -0.001, -0.001, -0.003, -0.005, -0.008, -0.010, -0.012, -0.014, -0.018, -0.023, -0.025, -0.065, -0.104, -0.180],
    [ -0.000, -0.001, -0.001, -0.002, -0.004, -0.006, -0.008, -0.010, -0.012, -0.015, -0.020, -0.031, -0.058, -0.093, -0.158],
    [ -0.000, -0.000, -0.001, -0.001, -0.003, -0.004, -0.006, -0.007, -0.009, -0.011, -0.015, -0.024, -0.046, -0.073, -0.122],
    [ -0.000, -0.000, -0.001, -0.001, -0.002, -0.003, -0.005, -0.006, -0.007, -0.009, -0.012, -0.020, -0.038, -0.060, -0.100]
]

SR1 = pd.DataFrame(data_SR1, columns=header, index=index)


# Change the name of the header and index
SR1 = SR1.rename_axis("Pr", axis=1)
SR1 = SR1.rename_axis("Tr", axis=0)
