import tkinter as tk
from tkinter import ttk

def start_timer():
    global remaining_time
    user_input = entry.get()
    if user_input.isdigit():
        remaining_time = int(user_input)
        update_timer()
    else:
        timer_label.config(text="Please enter a valid number")

def update_timer():
    global remaining_time
    if remaining_time > 0:
        timer_label.config(text=f"Time remaining: {remaining_time} seconds")
        remaining_time -= 1
        timer_label.after(1000, update_timer)
    else:
        timer_label.config(text="Time's up!")

# Create the main window
window = tk.Tk()
window.title("Timer App")

# Create and configure the Entry widget for user input
entry_label = ttk.Label(window, text="Set timer (seconds):")
entry_label.pack()
entry = ttk.Entry(window)
entry.pack()

# Create and configure the Start button
start_button = ttk.Button(window, text="Start Timer", command=start_timer)
start_button.pack()

# Create and configure the timer label
timer_label = ttk.Label(window, text="")
timer_label.pack()

# Initialize the remaining_time variable
remaining_time = 0

# Start the GUI event loop
window.mainloop()
