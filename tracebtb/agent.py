#!/usr/bin/env python

"""
    tracebtb agent
    Created Nov 2025
    Copyright (C) Damien Farrell

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 3
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
"""

import sys,os,subprocess,glob,re
import time, datetime
import platform
import pandas as pd
from langchain.tools import tool
import panel as pn
from panel.chat import ChatInterface
from . import tools

home = os.path.expanduser("~")
module_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(module_path,'data')
config_path = os.path.join(home, '.config','tracebtb')

@tool
def get_herd_movements(herd_id: str):
    """Fetches the last 12 months of cattle movements for a specific herd."""
    # Your code to query your database/dataframe
    return data

@tool
def calculate_genomic_distance(isolate_a: str, isolate_b: str):
    """Calculates the SNP distance between two M. bovis isolates."""
    # Your logic for SNP comparison
    return snp_dist

def test_chat_widget():
    """Test chat interface"""

    def even_or_odd(contents, user, instance: ChatInterface):
        if len(contents) % 2 == 0:
            return "Even number of characters."
        return "Odd number of characters."

    ch = ChatInterface(
        callback=even_or_odd,
        callback_user="Assistant",
        avatar=os.path.join(module_path, 'cow.png'),
        show_button_name=False
    )
    return ch

def ai_chat_widget(df):
    """
    Demonstrates how to use the ChatInterface widget to create an AI chatbot
    """

    from langchain_ollama import ChatOllama
    from langchain_experimental.agents.agent_toolkits import create_pandas_dataframe_agent

    #modelname = "deepseek-r1:14b"   
    modelname = 'gemma:7b'
    # We set temperature to 0 for precision in data analysis
    llm = ChatOllama(
        model=modelname,
        temperature=0,
        #stop=["<|file_separator|>"]
    )
    # Create the Agent
    agent = create_pandas_dataframe_agent(
        llm,
        df,
        verbose=True,
        allow_dangerous_code=True,
        #prefix=custom_prefix,
        agent_type="zero-shot-react-description"
    )

    async def callback(contents: str, user: str, instance: pn.chat.ChatInterface):
        response = agent.invoke(contents)
        message = response['output']
        yield message

    chat_interface = pn.chat.ChatInterface(callback=callback, callback_user="Assistant")
    chat_interface.send(
        "Send a message!", user="System", respond=False
    )
    return chat_interface