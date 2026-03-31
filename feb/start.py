#!/usr/bin/env python3
"""
Простой скрипт запуска с правильной настройкой импортов.
Использование: python start.py [аргументы]
"""

import sys
import os

# Добавляем feb в путь
feb_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, feb_dir)

# Запускаем main напрямую
exec(open(os.path.join(feb_dir, 'main.py')).read())
