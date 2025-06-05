#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import sys


if __name__ == "__main__":
    # Could have a metadata.out or something like that file output by the simulation code 
    # where I could just read this info in from. This would require more (unrequested)
    # file and string handling in C which makes me physically ill, so I won't do that.
    # For now, if needed, the user can pass these as commandline arguments. Could also 
    # consider calling animation python script from inside the C file.
    metadata = {
        "FILE": "",
        "GRID_HEIGHT": 500,
        "GRID_WIDTH": 500,
        "N_FRAMES": 500,
        "FPS": 20,
    }

    def parse_command_line_arguments():
        """
        Parse command line arguments

        This function checks the command line arguments passed after ./animate.py
        and matches them against FILE=<>, GRID_HEIGHT=<>, GRID_WIDTH=<>, N_FRAMES=<>, FPS=<>.
        If a match is found defining one of these variables, the default value stored in the
        metadata dictionary is overriden. The FILE parameter must always be passed.

        If the command line arguments are passed incorrectly, execution is stopped and the correct
        usage is printed to console.
        """

        if len(sys.argv) < 2:
            print("Usage:. /animate.py FILE = <relative_path_to_simulation_binary> GRID_HEIGHT=<> GRID_WIDTH=<> N_FRAMES=<> FPS=<>\n")
            print("Example: ./animate.py FILE = outputs/viscosity_investigation_13_05_2025_16h51m39s.bin")
            print("The parameters GRID_HEIGHT, GRID_WIDTH, N_FRAMES and desired animation FPS are optional parameters.\n")
            print("Defaults:\nGRID_HEIGHT=500\nGRID_WIDTH=500\nN_FRAMES=500\nFPS=20\n")
            print("If your simulation does not match these settings, you must pass them as commandline arguments.")
            sys.exit(1)

        for arg in sys.argv[1:]:
            if "=" in arg:
                key, value = arg.split("=", 1)
                if key in metadata.keys():
                    if key != "FILE":
                        metadata[key] = int(value)
                    else:
                        metadata[key] = value
                else:
                    print(f"Wrong key {key}. Must be one of {list(metadata.keys())}")
                    sys.exit(1)

    parse_command_line_arguments()

    try:
        data = np.fromfile(metadata["FILE"], dtype = np.float32).reshape(metadata["N_FRAMES"], metadata["GRID_HEIGHT"], metadata["GRID_WIDTH"])
        
        fig, ax = plt.subplots()
        im = ax.imshow(data[0], cmap="inferno", animated=True)
        fig.colorbar(im)
        
        def update(frame):
            im.set_array(data[frame])
            ax.set_title(f"Frame {frame}")
            return [im]

        print("Animating!\n")
        ani = animation.FuncAnimation(fig, update, frames=metadata["N_FRAMES"], blit=True)

        # Save as animated GIF (requires Pillow: sudo apt install python3-pillow)
        save_file = metadata["FILE"].removesuffix(".bin") + "_animation.gif"
        ani.save(save_file, writer='pillow', fps=metadata["FPS"])

    except ValueError:
        print("Binary file cannot be reshaped to (N_FRAMES, GRID_WIDTH, GRID_HEIGHT). One of these parameters is wrong.\n")
        print(f"The current values are:\nGRID_HEIGHT={metadata["GRID_HEIGHT"]}\nGRID_WIDTH={metadata["GRID_WIDTH"]}\nN_FRAMES={metadata["N_FRAMES"]}\nFPS={metadata["FPS"]}\n")
        print("If your simulation does not match these settings, you must pass them as commandline arguments.")
        sys.exit(1)

    except FileNotFoundError:
        print("Provided FILE is wrong.")
        sys.exit(1)

