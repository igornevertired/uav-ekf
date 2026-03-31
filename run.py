#!/usr/bin/env python3
"""
Скрипт запуска программы с правильной настройкой путей импорта.
"""

import sys
import os

# Добавляем родительскую директорию в PYTHONPATH
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

# Теперь импортируем главный модуль
if __name__ == "__main__":
    from feb.main import main
    main()
