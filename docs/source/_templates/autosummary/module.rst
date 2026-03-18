{{ fullname | escape | underline }}

.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :show-inheritance:

{% if classes %}
.. rubric:: Classes

.. autosummary::
   :toctree: .
   :template: autosummary/class.rst

{% for item in classes %}
   {{ item }}
{% endfor %}
{% endif %}

{% if functions %}
.. rubric:: Functions

.. autosummary::
   :toctree: .
   :template: autosummary/function.rst

{% for item in functions %}
   {{ item }}
{% endfor %}
{% endif %}
