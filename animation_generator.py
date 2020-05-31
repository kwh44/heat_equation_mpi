import imageio
import os

output_dir = "./results/"
filenames = os.popen("ls " + output_dir).read().split("\n")
with imageio.get_writer('./heatmap.gif', mode='I') as writer:
    for filename in filenames:
        if filename:
            try:
                image = imageio.imread(output_dir + filename)
                writer.append_data(image)
            except:
                print("Error: " + filename)
                continue
