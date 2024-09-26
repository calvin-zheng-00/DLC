import random
import pandas as pd

Dressing = {
    1:"Buttoning (break) Unbuttoning shirt (holes on dominant hand side)",
    2:"Buttoning (break) Unbuttoning studs (holes on dominant hand side)",
    3:"Opening/Closing zipper",
    4:"Opening/Closing velcro (Just the first one)",
    5:"Put on (break) take off sunglasses (orientated sideways to participant)"
}

Dining = {
    1:"Eat with spoon (face up, pointing up)",
    2:"Eat with fork (face down, pointing up)",
    3:"Saw with butter knife",
    4:"Pick up bowl (on tripod, back and forth is one set)",
    5:"Drink from cup with handle",
    6:"Drink from wine glass (stem)",
    7:"Drink from glass",
    8:"Drink from beer bottle",
    9:"Drink from mug",
    10:"Screw water bottle cap on (break) off",
    11:"Pick up plate (On tripod, back and forth is one set)",
    12:"Dip tea bag"
}

Cooking = {
    1:"Cut with big kitchen knife (blue)",
    2:"Grab pan (handle perpendicular to participant)",
    3:"Grab pan lid",
    4:"Use whisk",
    5:"Pick up strainer (handle horizontal participant)",
    6:"Pour jug",
    7:"Pour cereal",
    8:"Grab block at x with tongs, hold for one second, drop then return",
    9:"Use spatula (scraping off the table motion)",
    10:"Screw jar lid on (break) off"
}

Mobility = {
    1:"Turning door knob",
    2:"Turning door handle",
    3:"Using a joystick (one circle)",
    4:"Flick a switch (on and off in one motion)",
    5:"Hold umbrella upwards",
    6:"Pick up key and use key lock (ignoring reset)",
    7:"Twist dial lock (left then right)"
}

Bathing = {
    1:"Grabbing toilet paper roll",
    2:"Pulling tissue paper from the box and wiping the nose (drop at x)",
    3:"Wipe face with towel (laid out flat, ignoring reset)",
    4:"Brushing with a toothbrush",
    5:"Grab toothpaste",
    6:"Twist toothpaste cap on (break) off",
    7:"Use hair comb",
    8:"Use hair brush",
    9:"Use hair dryer (need to turn it on with button press)"
}

Cleaning = {
    1:"Pick up sponge from one location then wipe with it in another before returning it",
    2:"Use cleaning spray",
    3:"Pick up brush from one location then wipe with it in another before returning it",
    4:"Use iron"
}

Writing = {
    1:"Click pen on and off",
    2:"Writing with a pencil (Hello (draw/write from left to right, each activity on unique line))",
    3:"Writing with a pencil (Circle the size of a 50 cent coin (draw/write from left to right, each activity on unique line))",
    4:"Using an eraser",
    5:"Pick up paper",
    6:"Pick up a ruler",
    7:"Open and flip page of book",
    8:"Use scissors",
    9:"Pick up card",
    10:"Use stapler (Keep on table, press down)",
    11:"pick up binder clip and open",
    12:"Pick up USB",
    13:"Use tape dispenser (place tape at x)"
}

Digital = {
    1:"Typing on a keyboard (Hello)",
    2:"Trace a circle with a mouse",
    3:"Mouse click and scroll (double left click, double right click, scroll up, scroll down)",
    4:"Type on a calculator (2+2=) (Don't pick up)",
    5:"Pick up mobile phone to call",
    6:"Type on phone (Hello) (pick up)",
    7:"Pick up and press remote (one button once, aim at camera)",
    8:"Replace batteries in remote (place in (break) take out)"
}

Physical = {
    1:"Use screwdriver (do three screws on imaginary screw at x)",
    2:"Use power Drill (on imaginary screw at x)",
    3:"Use hammer (on towel placed at x, three hits)",
    4:"Pick up tape",
    5:"Connecting (break) disconnecting power plug (Use disconnected head pointing right, board with earth pointing down)",
    6:"Use wrench (do three screws on imaginary nut at x)",
    7:"Pick up small bucket",
    8:"Pick up screw",
}

General = {
    1:"Peg board (break, capture reset just in case)",
    2:"Pick up block",
    3:"Stacking coins (break, capture reset just in case) (reset one coin at a time)",
    4:"Pick up ball"
}

Activities = {
    1:Dressing,
    2:Dining,
    3:Cooking,
    4:Mobility,
    5:Bathing,
    6:Cleaning,
    7:Writing,
    8:Digital,
    9:Physical,
    10:General
}

Activities_name = {
    1:"Dressing",
    2:"Dining",
    3:"Cooking",
    4:"Mobility",
    5:"Bathing",
    6:"Cleaning",
    7:"Writing",
    8:"Digital",
    9:"Physical",
    10:"General"
}

f = open("order.txt", "w")
inst_list = []

category_order = list(range(1,len(Activities)+1))
random.shuffle(category_order)
for category in category_order:
    f.write(Activities_name[category])
    f.write(':\n')
    activity_order = list(range(1,len(Activities[category])+1))
    random.shuffle(activity_order)
    for activity in activity_order:
        f.write('\t')
        f.write(Activities[category][activity])
        f.write('\n')
        inst_str = Activities_name[category][:3]
        inst_str = inst_str.upper()
        inst_str = "Instructions_"+inst_str+str(activity)
        inst_list.append(inst_str)
        

f.close()
inst_df = pd.DataFrame(inst_list, columns=['videos'])
inst_df.to_csv("C:/Users/czhe0008/Documents/DLCprojects/openpose/data_conversion/activity_order/instructions.csv", index = False)