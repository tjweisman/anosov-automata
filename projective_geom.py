import numpy as np

class ProjectivizationException:
    pass

def affine_coords(points, chart_index = None):

    apoints = np.array(points)
    
    _chart_index = chart_index
    if chart_index is None:
        _chart_index = np.argmax(np.amin(np.absolute(points), axis=0))

    if (apoints[:,_chart_index] == 0).any():
        if chart_index is not None:
            raise ProjectivizationException(
                "points don't lie in the specified affine chart"
            )
        else:
            raise ProjectivizationException(
                "points don't lie in any standard affine chart"
            )
    
    normalized = (apoints / apoints[:,_chart_index][:,None]).T
    
    affine = np.concatenate(
        [normalized[:_chart_index], normalized[_chart_index + 1:]]).T


    if chart_index is None:
        return (affine, _chart_index)
    else:
        return affine

def lift(points, chart_index):
    return np.insert(points, chart_index, np.ones(points.shape[0]), axis=1)

def change_chart(points, chart_index, new_index):
    return affine_coords(lift(points, chart_index), new_index)

def sphere_lift(points):
    return lift(*affine_coords(points))
