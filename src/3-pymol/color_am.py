from pymol import cmd

def coloram(**kwargs):
    """
    AUTHOR
    Tamas Hegedus @ hegelab.org

    DESCRIPTION
    Colors structures according to b-factors with AlphaMissense colors.

    USAGE
    coloram sel=<selection>

    PARAMETERS
    sel (string)
    The name of the selection/object to color with AlphaMissense colors. If not provided, an error message will be displayed.
    """

    # Default Selection
    selection = str(kwargs.get("sel"))
    
    if not selection or selection.lower() == "none":
        print("Error: No sel argument provided")
        return

    # Define colors
    cmd.set_color("benign", [161, 176, 222])
    cmd.set_color("ambiguous", [230, 230, 230])
    cmd.set_color("patho0", [245, 125, 129])
    cmd.set_color("patho1", [165, 13, 18])

    # Apply colors based on b-factors
    cmd.color("benign", f"({selection}) and b < 0.34")
    cmd.color("ambiguous", f"({selection}) and (b = 0.34 or b > 0.34)")
    cmd.color("patho0", f"({selection}) and (b = 0.564 or (b > 0.564 and b < 0.78))")
    cmd.color("patho1", f"({selection}) and (b = 0.78 or b > 0.78)")

# Extend PyMOL with the function
cmd.extend("coloram", coloram)
cmd.auto_arg[0]["coloram"] = [cmd.object_sc, "object", ""]
