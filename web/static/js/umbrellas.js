$(function(){

    // Replica model
    window.Replica = Backbone.Model.extend({
        defaults: {
            "name":  "new_replica",
            "parameters": { coordinate: "", coordinates: "", force: "" },
        },
        
        addParameter: function(k,v) {
            this.parameters[k] = v;
        },
        
        removeParameter: function(k) {
        },
        
        getParameter: function(k) {
            return this.parameters[k]
        },
        
    });

    // ReplicaSet
    window.ReplicaSet = Backbone.Collection.extend({
        model: Replica,
        url: '/replicas',
    });

    window.Replicas = new ReplicaSet;
    
    // Replica view
    window.ReplicaView = Backbone.View.extend({
        tagName: "li",
        template: _.template($('#replica-template').html()),
        
        events: {
        },
        
        initialize: function() {
            _.bindAll(this, 'render');
            // Bind the view to the model (one-to-one)
            this.model.bind('change', this.render);
            this.model.view = this;
        },
        
        // Render the contents of a single replica.
        render: function() {
            $(this.el).html(this.template(this.model.toJSON()));
            return this;
        },
        
    });
    
    // The Application

    // Our overall **AppView** is the top-level piece of UI.
    window.AppView = Backbone.View.extend({
        el: $("#umbrellas-app"),

        events: {
        },

        initialize: function() {
            _.bindAll(this, 'addOne', 'addAll', 'render');
            Replicas.bind('add',     this.addOne);
            Replicas.bind('refresh', this.addAll);
            Replicas.bind('all',     this.render);
        },

        render: function() {
            //this.$('#ensemble-stats').html(this.chartTemplate({ }));
        },

        // Add a single replica to the list by creating a view for it, and
        // appending its element to the `<ul>`.
        addOne: function(replica) {
            var view = new ReplicaView({model: replica});
            this.$("#replica-list").append(view.render().el);
        },
        
        addAll: function() {
            Replicas.each(this.addOne);
        }

    });

    // Finally, we kick things off by creating the **App**.
    window.App = new AppView;

});