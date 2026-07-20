# Уровни многословности консольного вывода

## Уровни вывода (Verbosity Levels)

Консольный вывод симуляции контролируется параметром `verbosity_level`:

| Level | Описание | Используется для |
|-------|----------|------------------|
| **0** | 🔇 Полный silence | Батч-прогоны без console spam; фокус на диагностике в файлы |
| **1** | 🟢 Step markers (default) | Нормальные прогоны: видно прогресс, но не шум о каждой операции |
| **2** | 🔍 Verbose: все таймеры | Отладка, поиск узких мест, детальное профилирование |

## Примеры использования

### В Lua (перед `require("sim_core").run()`)

```lua
-- Тихий прогон
pic_config.set_verbosity_level(0)
require("sim_core").run()

-- Нормальный прогон (default)
pic_config.set_verbosity_level(1)
require("sim_core").run()

-- Детальный вывод
pic_config.set_verbosity_level(2)
require("sim_core").run()
```

## Примеры вывода

**Level 0 (silent):**
```
OpenPIC: Simulation started ...
   0 sec: scatter NP and UP at t=0
  12 sec: Step 0: Save initial state
  ...
  [диагностика идёт в файлы, консоль чиста]
```

**Level 1 (default, step markers only):**
```
Step 1 of 1442
Step 2 of 1442
Step 3 of 1442
...
Step 100 of 1442
```

**Level 2 (verbose, all timers):**
```
Step 1 of 1442
 363 sec: calc_magnetic_field_half_time
 364 sec: calc_electrons_velocity
 364 sec: update_Spitzer_coefficients
 365 sec: HeatSolver::solve
 365 sec: calc_electric_field
 366 sec: move_particles_half_time
 ...
Step 2 of 1442
```

## Рекомендации

| Сценарий | Level |
|----------|-------|
| Production batch runs (кластер, много кейсов подряд) | **0** |
| Interactive development, живой мониторинг | **1** (default) |
| Отладка, поиск bottleneck | **2** |
| Верификационные/benchmark прогоны | **1** или **0** |

## Технические детали

- **Level 0** подавляет ВСЕ `print_tm()` вызовы
- **Level 1** выводит только `Step X of Y` на каждом шагу (без промежуточных таймеров)
- **Level 2** выводит всё: шаги + промежуточные таймеры операций (calc, boundary conditions и т.д.)
- Каждый `print_tm()` вызов дополнительно проверяет что это процесс 0 в MPI режиме

## Изменение на ходу

Параметр можно изменить в любой момент **перед** `require("sim_core").run()`:

```lua
-- Начальное значение
pic_config.set_verbosity_level(2)  -- verbose

-- Потом изменить
pic_config.set_verbosity_level(0)  -- silent

-- Запустить
require("sim_core").run()
```
